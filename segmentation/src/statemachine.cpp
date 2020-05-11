
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>

const int HUNT=1;
const int SURFACE=2;
const int DIVE=3;
const int CLIMB=4;
const int DRIFT=5;

struct Segment {
    int start, end, length;
    int label;
    int state;
    Segment() {}
    Segment(double l[5]) {
        start=l[0];
        end=l[1];
        length=l[2];
        label=l[3];
        state=l[4];
    }
};

int main(int argc, char * argv[]) {
    assert(argc>1);
    FILE *fp=fopen(argv[1],"r");
    assert(fp);
    std::vector<Segment> seg;
    while (!feof(fp)) {
        double l[5];
        char buffer[1024] = {0,};
        if (fgets(buffer,1024,fp)==NULL) {
            continue;
        }
        if (sscanf(buffer," %le %le %le %le %le ",l+0,l+1,l+2,l+3,l+4) != 5) {
            continue;
        }
        seg.push_back(Segment(l));
    }
    fclose(fp);
    int state = -1;

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
                    state = 2;
                } else if (seg[i].label!=SURFACE) {
                    state=1;
                }
                break;
            case 1:
                if (seg[i].label==DIVE) {
                    state = 2;
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
                if (seg[N-i].label==DIVE) {
                    state = 2;
                } else if (seg[N-i].label!=SURFACE) {
                    state=1;
                }
                break;
            case 1:
                if (seg[N-i].label==DIVE) {
                    state = 2;
                }
                break;
            case 2:
                if (seg[N-i].label!=DIVE) {
                    state = -1;
                }
                break;
            default:
                break;
        }
        seg[N-i].state = state;
        if (state==1) {
            seg[N-i].label = DIVE;
        } else if ((state == -1) && (seg[N-i].label==DIVE)) {
            seg[N-i].label = HUNT;
        }
    }

    
    printf("Output\n");
    fp=fopen(argv[1],"w");
    for (size_t i=0;i<seg.size();i++) {
        fprintf(fp,"%d %d %d %d %d\n",
                seg[i].start,seg[i].end,seg[i].length,
                seg[i].label,seg[i].state);
    }
    fclose(fp);

    return 0;
}

