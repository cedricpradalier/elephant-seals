/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#include <math.h>
#include <signal.h>

#include "cgnuplot/CGnuplot.h"
#include "cgnuplot/CGnuplotSurfaceGrapher.h"

using namespace cgnuplot;
bool signal_received = false;

void sighdl(int n)
{
	signal_received = true;
}

const int size = 23;

void fillmatrix(vector<double>& m, int v)
{
    m.clear();
    int k=0;
    for (int i=0; i<size; ++i)
    {
	for (int j=0; j<size; ++j)
	{
	    m.push_back(cos((double) (i*i + j*j)*100. / (size*size*(v+1.))));

	    ++k;
	}
    }
}

void fillmatrix2(vector<double>& m, int v)
{
    m.clear();
    int k=0;
    for (int i=0; i<size; ++i)
    {
	for (int j=0; j<size; ++j)
	{
	    m.push_back(sin((double) (cos((double) i*i)*j*j)*100. / (size*size*(v+1.))));

	    ++k;
	}
    }
}

int main(int argc, char * argv[])
{
	int i;
	// Le signal est défini apres le fork, pour qu'il ne concerne que ce
	// process, et pas gnuplot.
	signal(SIGINT,sighdl);

	
	// Crée une nouvelle instance de gnuplot (fork...)
	// Binaire par defaut et options par defaut.
	CGnuplot gpl;
	CGnuplotSurfaceGrapher surf(gpl);
	//gpl.setLogging(true);
	//gpl.setEcho(true);

	surf.init_gnuplot(50, 50, true);

	vector<double> m;

	for (i=1;i<200;i++)
	{
		if (signal_received) break;
		// Affichage de la fonction. La fonction rend la main quand gnuplot
		// rend la main...

		fillmatrix(m, 5+abs(50 - i%100));
		surf.draw(m, size, -1, 1, "le joli test de surface!");
		fillmatrix2(m, 5+abs(50 - i%100));
		surf.draw(m, size, -1, 2, "rooooh, une autre fenetre en parallele!");
	}
	
	return 0;
}
