
#include <pthread.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <complex.h>
#include <fftw3.h>


const unsigned int MAXFFT = 5000;
const unsigned int NTHREADS = 32;
const int HUNT=1;
const int SURFACE=2;
const int DIVE=3;
const int CLIMB=4;
const int DRIFT=5;

struct DataPoint {
    double t;
    float A[3];
    float P, dP;
    DataPoint() {}
    DataPoint(double l[5]) {
        t=l[0];
        A[0]=l[1];
        A[1]=l[2];
        A[2]=l[3];
        P=l[4];
        dP=0;
    }
    float An() const {
        return sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    }
};

struct FFTPlan {
    size_t N;
    double *in;
    fftw_complex *out;
    fftw_plan p;

    FFTPlan(size_t N) : N(N) {
        in = (double*) fftw_malloc(sizeof(double) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);
    }

    ~FFTPlan() {
        fftw_destroy_plan(p);
        fftw_free(in); fftw_free(out);
    }

    void clear() {
        for (size_t i=0;i<N;i++) {
            out[i][0]=out[i][1]=in[i]=0;
        }
    }

    void execute(const std::vector<DataPoint> & data) {
        clear();
        double mean = 0;
        for (size_t i=0;i<std::min(N,data.size());i++) {
            in[i] = data[i].A[1];
            mean += in[i];
        }
        mean /= data.size();
        for (size_t i=0;i<std::min(N,data.size());i++) {
            in[i] -= mean;
        }
        fftw_execute(p); 
    }
};

    

struct Segment {
    size_t start, end, length, diveid;
    int label;
    int state;
    std::vector<DataPoint> data;
    cv::Mat_<float> fft;
    Segment() {}
    Segment(double l[5]) {
        start=l[0];
        end=l[1];
        length=l[2];
        label=l[3];
        state=l[4];
	diveid=0;
    }
    void updateData(const std::vector<DataPoint> & all) {
        for (size_t i=start;i<=end;i++) {
            data.push_back(all[i]);
        }
	assert(data[1].t > data[0].t);
    }

    void computeFFT(FFTPlan & P) {
        
        P.execute(data);
        fft = cv::Mat_<float>(P.N,5,0.0);
        float dt = (data[1].t-data[0].t);
	assert((dt > 0) && (dt < 1));
	dt *= P.N;
	// printf("Sampling %.2f freq step %.3f\n",(data[1].t-data[0].t),dt);
        for (size_t i=0;i<P.N;i++) {
            fft(i,0) = i/dt;
            fft(i,1) = P.in[i];
            fft(i,2) = P.out[i][0];
            fft(i,3) = P.out[i][1];
            fft(i,4) = hypot(P.out[i][0],P.out[i][1]);
#if 0
            printf("%d %.2e %.2e %.2e %.2e %.2e\n",
                    int(i),fft(i,0),fft(i,1),fft(i,2),
                    fft(i,3),fft(i,4));
#endif
        }
        //printf("\n");
    }

    double depthMax() const {
        double dmax=0;
        for(size_t i=0;i<data.size();i++) {
            dmax=std::min<double>(data[i].P,dmax);
        }
        return dmax;
    }


};

double estimateDiveVelocity(const std::vector<DataPoint> & data, 
        unsigned int i0, unsigned int i1) {
    cv::Mat_<float> A(i1-i0+1,2,1.0f);
    cv::Mat_<float> B(i1-i0+1,1,0.0f);
    for (unsigned int i=i0;i<=i1;i++) {
        A(i-i0,0) = data[i].t;
        B(i-i0,0) = data[i].P;
    }
    cv::Mat_<float> X;
    cv::solve(A, B, X , cv::DECOMP_QR);
    if (isnan(X(0,0))) return 0;
    return X(0,0);
}

struct ThreadContext {
	pthread_mutex_t mutex;
	size_t index;
	std::vector<DataPoint> & data;

	ThreadContext(std::vector<DataPoint> & data) : data(data) {
		index = 0;
		pthread_mutex_init(&mutex,NULL);
	}
	~ThreadContext() {
		pthread_mutex_destroy(&mutex);
	}

	void run() {
		while (1) {
			size_t i = data.size();
			pthread_mutex_lock(&mutex);
			i = index;
			index += 1;
			pthread_mutex_unlock(&mutex);
			if (i >= data.size()) {
				break;
			}
			data[i].dP = estimateDiveVelocity(data,std::max<int>(0,int(i)-10),
					std::min<int>(data.size()-1,int(i)+10));
		}
	}
};

void * thread_fun(void * arg) {
	ThreadContext *context = (ThreadContext*)arg;
	context->run();
	return NULL;
}


int main(int argc, char * argv[]) {
    assert(argc>2);
    std::map<std::string,int> label_name;
    label_name["SURFACE"]=SURFACE;
    label_name["DIVE"]=DIVE;
    label_name["HUNT"]=HUNT;
    label_name["CLIMB"]=CLIMB;
    label_name["DRIFT"]=DRIFT;
    std::vector<DataPoint> data;
    std::vector<unsigned char> labels;
    std::vector<Segment> seg;
    FILE *fp=fopen(argv[1],"r");
    assert(fp);
    while (!feof(fp)) {
        int v;
        char buffer[1024] = {0,};
        if (fgets(buffer,1024,fp)==NULL) {
            continue;
        }
        if (sscanf(buffer," %d ",&v) != 1) {
            continue;
        }
        labels.push_back(v);
    }
    fclose(fp);
    printf("Loaded %d labels\n",int(labels.size()));

    fp=fopen(argv[2],"r");
    assert(fp);
    while (!feof(fp)) {
        char buffer[1024] = {0,};
        if (fgets(buffer,1024,fp)==NULL) {
            continue;
        }
        double l[5];
#if 0
        if (sscanf(buffer," %e , %e , %e , %e , %e ",l+0,l+1,l+2,l+3,l+4) != 5) {
            printf("Rejecting '%s'\n",buffer);
            continue;
        }
#else
	long int v[5]={0,0,0,0,0};
        if (sscanf(buffer," %ld %ld %ld %ld %ld ",v+0,v+1,v+2,v+3,v+4) != 5) {
            printf("Rejecting '%s'\n",buffer);
            continue;
        }
	for (int i=0;i<5;i++) {
		l[i] = double(v[i]/100) + (v[i]%100)/100.0;
	}
	l[4] = -l[4];
	if (data.size()>0) {
		assert(l[0] >= data[data.size()-1].t);
	}
#endif
        data.push_back(DataPoint(l));
    }
    fclose(fp);
    printf("Loaded %d datapoints\n",int(data.size()));
    assert(labels.size() == data.size());

    mkdir("output",0755);

#if 1
    pthread_t tid[NTHREADS];
    ThreadContext context(data);
    for (unsigned int i=0;i<NTHREADS;i++) {
	    pthread_create(tid+i,NULL,thread_fun,&context);
    }
    for (unsigned int i=0;i<NTHREADS;i++) {
	    pthread_join(tid[i],NULL);
    }
#else
    for (size_t i=0;i<data.size();i++) {
        data[i].dP = estimateDiveVelocity(data,std::max<int>(0,int(i)-10),
                std::min<int>(data.size()-1,int(i)+10));
    }
#endif
    printf("Computed velocity\n");

    Segment segment;
    segment.start=0;
    for (size_t i=1;i<labels.size();i++) {
        if (labels[i]!=labels[i-1]) {
            segment.end=i-1;
            segment.length = segment.end-segment.start+1;
            segment.label = labels[i-1];
            segment.state = 0;
            seg.push_back(segment);
            segment.start=i;
        }
    }
    segment.end = labels.size()-1;
    segment.length = segment.end-segment.start+1;
    segment.label = labels[labels.size()-1];
    seg.push_back(segment);


    std::vector<Segment> segf;
    for (size_t i=0;i<seg.size();i++) {
        Segment s = seg[i];
        if (s.length > 50) {
            if (segf.size()==0) {
                s.start = 0;
            }
            segf.push_back(s);
        } else if (segf.size()>0) {
            segf[segf.size()-1].end = s.end;
        }
    }
    segf[segf.size()-1].end = seg[seg.size()-1].end;
    seg=segf;
    segf.clear();
    int state = -1;
    
#if 1
    printf("Backward pass\n");
    state = -1;
    for (size_t i=0;i<seg.size();i++) {
        size_t N=seg.size() - 1;
        // Transitions
        switch (state) {
            case -1:
                if (seg[N-i].label==SURFACE) {
                    state = 0;
                }
                break;
            case 0:
                if (seg[N-i].label==CLIMB) {
                    state = 2;
                } else if (seg[N-i].label!=SURFACE) {
                    state=1;
                }
                break;
            case 1:
                switch (seg[N-i].label) {
                    case DIVE:
                        state = -1;
                        break;
                    case CLIMB:
                        state = 2;
                        break;
                    case SURFACE:
                        state = 0;
                        break;
                    default:
                        break;
                }
                break;
            case 2:
                if (seg[N-i].label!=CLIMB) {
                    state = -1;
                }
                break;
            default:
                break;
        }
        seg[N-i].state = state;
        if (state==1) {
            seg[N-i].label = CLIMB;
        } else if ((state == -1) && (seg[N-i].label==CLIMB)) {
            seg[N-i].label = HUNT;
        }
    }
#endif


#if 1
    // forward pass
    printf("Forward pass\n");
    state = -1;
    for (size_t i=0;i<seg.size();i++) {
        // Transitions
        switch (state) {
            case -1:
                if (seg[i].label==SURFACE) {
                    state = 0;
                }
                break;
            case 0:
                if (seg[i].label==DIVE) {
                    state=2;
                } else if (seg[i].label!=SURFACE) {
                    state=1;
                }
                break;
            case 1:
                switch (seg[i].label) {
                    case DIVE:
                        state = 2;
                        break;
                    case CLIMB:
                        state = -1;
                        break;
                    case SURFACE:
                        state = 0;
                        break;
                    default:
                        break;
                }
                break;
            case 2:
                if (seg[i].label!=DIVE) {
                    state = -1;
                }
                break;
            default:
                break;
        }
        seg[i].state = state;
        if (state==1) {
            seg[i].label = DIVE;
        } else if ((state == -1) && (seg[i].label==DIVE)) {
            seg[i].label = HUNT;
        }
    }
#endif


    // Regenerating labels
    printf("Regenerating labels\n");
    std::vector<unsigned char> labels_orig(labels);
    labels.clear();
    for (size_t i=0;i<seg.size();i++) {
        seg[i].length = seg[i].end - seg[i].start + 1;
        for (size_t j=0;j<(unsigned)seg[i].length;j++) {
            labels.push_back(seg[i].label);
        }
    }
    seg.clear();
    size_t dive_count = 0;
    segment.start=0;
    for (size_t i=1;i<labels.size();i++) {
        if (labels[i]!=labels[i-1]) {
            segment.end=i-1;
            segment.length = segment.end-segment.start+1;
            segment.label = labels[i-1];
            segment.state = 0;
	    segment.diveid = dive_count;
            seg.push_back(segment);
            segment.start=i;
	    if (labels[i]==SURFACE) {
		    dive_count += 1;
	    }

        }
    }
    segment.end = labels.size()-1;
    segment.length = segment.end-segment.start+1;
    segment.label = labels[labels.size()-1];
    segment.diveid = dive_count;
    seg.push_back(segment);
    
    typedef std::vector<size_t> IndexVector;
    std::map<int,IndexVector> DataSegments;
    for (size_t i=0;i<seg.size();i++) {
        seg[i].updateData(data);
        DataSegments[seg[i].label].push_back(i);
        
    }

#if 1
    printf("Generating FFTs: %d climbs\n",int(DataSegments[CLIMB].size()));
    size_t max_climb_size = 0;
    for (size_t i=0;i<DataSegments[CLIMB].size();i++) {
        size_t j = DataSegments[CLIMB][i];
        max_climb_size = std::max(max_climb_size,seg[j].data.size());
    }
    max_climb_size = std::min<size_t>(max_climb_size,MAXFFT);
    {
	    FFTPlan plan(max_climb_size);
	    // cv::Mat_<float> spectroinput(plan.N,DataSegments[CLIMB].size(),0.0);
	    cv::Mat_<float> spectrogram(plan.N,DataSegments[CLIMB].size());
	    for (size_t i=0;i<DataSegments[CLIMB].size();i++) {
		    size_t j = DataSegments[CLIMB][i];
		    seg[j].computeFFT(plan);
		    if ((i%100)==0) {
			    printf("+");fflush(stdout);
		    }
		    for (int k=0;k<spectrogram.rows;k++) {
			    spectrogram(k,i) = seg[j].fft(k,4);
		    }
		    // for (int k=0;k<spectroinput.rows;k++) {
		    //     spectroinput(k,i) = seg[j].fft(k,1);
		    // }
	    }
	    printf("\n");
	    fp = fopen("output/spectro_climb.dat","w");
	    fprintf(fp,"0 ");
	    for (int i=0;i<spectrogram.cols;i++) {
		    size_t j = DataSegments[CLIMB][i];
		    fprintf(fp,"%.2f ",seg[j].depthMax());
	    }
	    fprintf(fp,"\n");
	    fprintf(fp,"0 ");
	    for (int i=0;i<spectrogram.cols;i++) {
		    size_t j = DataSegments[CLIMB][i];
		    fprintf(fp,"%.2f ",seg[j].data[0].t);
	    }
	    fprintf(fp,"\n");
	    // fprintf(fp,"0 ");
	    // for (int i=0;i<spectrogram.cols;i++) {
	    //         size_t j = DataSegments[CLIMB][i];
	    //         fprintf(fp,"%lu ",seg[j].diveid);
	    // }
	    // fprintf(fp,"\n");
	    for (int i=0;i<spectrogram.rows;i++) {
		    size_t j = DataSegments[CLIMB][0];
		    double f = seg[j].fft(i,0);
		    if (f>2) break;
		    fprintf(fp,"%.3e ",f);
		    for (int j=0;j<spectrogram.cols;j++) {
			    fprintf(fp,"%.3e ", spectrogram(i,j));
		    }
		    fprintf(fp,"\n");
	    }
	    fclose(fp);
    }
#endif

#if 1
    printf("Generating FFTs: %d dive\n",int(DataSegments[DIVE].size()));
    size_t max_dive_size = 0;
    for (size_t i=0;i<DataSegments[DIVE].size();i++) {
        size_t j = DataSegments[DIVE][i];
        max_dive_size = std::max(max_dive_size,seg[j].data.size());
    }
    max_dive_size = std::min<size_t>(max_dive_size,MAXFFT);
    {
	    FFTPlan plan(max_dive_size);
	    cv::Mat_<float> spectrogram(plan.N,DataSegments[DIVE].size());
	    for (size_t i=0;i<DataSegments[DIVE].size();i++) {
		    size_t j = DataSegments[DIVE][i];
		    seg[j].computeFFT(plan);
		    if ((i%100)==0) {
			    printf("+");fflush(stdout);
		    }
		    for (int k=0;k<spectrogram.rows;k++) {
			    spectrogram(k,i) = seg[j].fft(k,4);
		    }
		    // for (int k=0;k<spectroinput.rows;k++) {
		    //     spectroinput(k,i) = seg[j].fft(k,1);
		    // }
	    }
	    printf("\n");
	    fp = fopen("output/spectro_dive.dat","w");
	    fprintf(fp,"0 ");
	    for (int i=0;i<spectrogram.cols;i++) {
		    size_t j = DataSegments[DIVE][i];
		    fprintf(fp,"%.2f ",seg[j].depthMax());
	    }
	    fprintf(fp,"\n");
	    fprintf(fp,"0 ");
	    for (int i=0;i<spectrogram.cols;i++) {
		    size_t j = DataSegments[DIVE][i];
		    fprintf(fp,"%.2f ",seg[j].data[0].t);
	    }
	    fprintf(fp,"\n");
	    // fprintf(fp,"0 ");
	    // for (int i=0;i<spectrogram.cols;i++) {
	    //         size_t j = DataSegments[DIVE][i];
	    //         fprintf(fp,"%lu ",seg[j].diveid);
	    // }
	    // fprintf(fp,"\n");
	    for (int i=0;i<spectrogram.rows;i++) {
		    size_t j = DataSegments[DIVE][0];
		    double f = seg[j].fft(i,0);
		    if (f>2) break;
		    fprintf(fp,"%.3e ",f);
		    for (int j=0;j<spectrogram.cols;j++) {
			    fprintf(fp,"%.3e ", spectrogram(i,j));
		    }
		    fprintf(fp,"\n");
	    }
	    fclose(fp);
    }
#endif

    
    printf("Output\n");
    


    fp=fopen("output/seg.csv","w");
    for (size_t i=0;i<seg.size();i++) {
        fprintf(fp,"%lu %lu %lu %lu %d %d %.1f %.2f %.1f %.2f %.2f\n",
                seg[i].diveid,seg[i].start+1,seg[i].end+1,seg[i].length,
                seg[i].label,seg[i].state,
                data[seg[i].start].t, data[seg[i].start].P,
                data[seg[i].end].t, data[seg[i].end].P,seg[i].depthMax());
    }
    fclose(fp);

    fp=fopen("output/labels.csv","w");
    for (size_t i=0;i<labels.size();i++) {
        fprintf(fp,"%d\n",labels[i]);
    }
    fclose(fp);

#if 0
    size_t i0=0,k=0;
    for (size_t i=0;i<seg.size();i++) {
        if ((i<seg.size()-1) && (seg[i].label != SURFACE)) continue;
        if ((i==seg.size()-1) || ((seg[i].end - i0) > 24*60*60*5)) {
            char fname[1024];
            sprintf(fname,"output/data%03d.csv",int(k));
            size_t i1=(i==seg.size()-1)?seg[i].end:((seg[i].start+seg[i].end)/2);
            fp=fopen(fname,"w");
            for (size_t i=i0;i<i1;i++) {
                fprintf(fp,"%.2f %.2f %.2f %.2f %.2f %.2f %d\n",
                        data[i].t,data[i].A[0],data[i].A[1],data[i].A[2],
                        data[i].An(), data[i].P, labels[i]);
            }
            fclose(fp);
            i0=i1;
            k+=1;
        }
    }
#endif

#if 0
    for (std::map<std::string,int>::const_iterator it=label_name.begin();
            it!=label_name.end();it++) {
        std::string dir=("output/"+it->first);
        mkdir(dir.c_str(),0755);
        for (size_t i=0;i<seg.size();i++) {
            if (seg[i].label != it->second) continue;
            char fname[1024];
            sprintf(fname,"%s/data%06d.csv",dir.c_str(),int(i));
            fp=fopen(fname,"w");
            for (size_t j=seg[i].start;j<=seg[i].end;j++) {
                fprintf(fp,"%.2f %.2f %.2f %.2f %.2f %.2f %.2f %d\n",
                        data[j].t,data[j].A[0],data[j].A[1],data[j].A[2],
                        data[j].An(), data[j].P, data[j].dP, labels[j]);
            }
            fclose(fp);
        }
    }
#endif

#if 1
    fp=fopen("output/dives.csv","w");
    for (size_t i=0;i<seg.size();i++) {
        if (seg[i].label != DIVE) continue;
        // Remove shallow dive. They are not relevant here
        if (seg[i].depthMax()>-150) continue;
        for (size_t j=seg[i].start;j<=seg[i].end;j++) {
            if (isnan(data[j].A[2])) continue;
            if (data[j].A[2]<0) continue;
            if (data[j].A[2]>10) continue;
            if (data[j].dP>0) continue;
            fprintf(fp,"%lu %.2f %.3f %.3f\n",seg[i].diveid,data[seg[i].start].t,
                    data[j].A[2],data[j].dP);
        }
    }
    fclose(fp);
#endif

#if 1
    fp=fopen("output/divesum.csv","w");
    double dacc = 3*3600; // accumulation duration
    double tacc = -1;
    unsigned int hullcount = 0;
    std::vector<cv::Point2f> accumulator;
    std::vector<cv::Point2f> hull;
    accumulator.push_back(cv::Point2f(0,1));
    accumulator.push_back(cv::Point2f(10,1));
    for (size_t i=0;i<seg.size();i++) {
        float tseg = data[seg[i].start].t;
        if (seg[i].label != DIVE) continue;
        // Remove shallow dive. They are not relevant here
        if (seg[i].depthMax()>-150) continue;
        if (tacc < 0) {
            tacc = tseg;
        }
        if (tseg - tacc > dacc) {
            hull.clear();
            hullcount += 1;
            cv::convexHull(accumulator,hull);
            for (size_t j=0;j<hull.size();j++) {
                if (hull[j].y>0) continue;
                fprintf(fp,"%.2f %.3f %.3f\n",tacc,hull[j].x,hull[j].y);
            }
            fprintf(fp,"\n");
            tacc = tseg;
            accumulator.clear();
            accumulator.push_back(cv::Point2f(0,1));
            accumulator.push_back(cv::Point2f(10,1));
        }
        for (size_t j=seg[i].start;j<=seg[i].end;j++) {
            if (isnan(data[j].A[2])) continue;
            if (data[j].A[2]<0) continue;
            if (data[j].A[2]>10) continue;
            if (data[j].dP>0) continue;
            accumulator.push_back(cv::Point2f(data[j].A[2],data[j].dP));
        }
    }
    if (!accumulator.empty()) {
        hull.clear();
        hullcount += 1;
        cv::convexHull(accumulator,hull);
        for (size_t j=0;j<hull.size();j++) {
            if (hull[j].y>0) continue;
            fprintf(fp,"%.2f %.3f %.3f\n",tacc,hull[j].x,hull[j].y);
        }
    }
    fclose(fp);
    printf("Done %d hulls\n",int(hullcount));
#endif

    return 0;
}

