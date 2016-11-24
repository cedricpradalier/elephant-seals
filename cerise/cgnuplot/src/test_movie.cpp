/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#include <signal.h>
#include <unistd.h>

#include "cgnuplot/CGnuplot.h"

using namespace cgnuplot;
int signal_received = 0;
CGnuplot * gpl;

void sighdl(int n)
{
	if (signal_received>1) return;
	if (signal_received>=1) {gpl->interrupt();}
	signal_received += 1;
	
}

int main(int argc, char * argv[])
{
	int i;
	
	// Crée une nouvelle instance de gnuplot (fork...)
	// Binaire par defaut et options par defaut.
	gpl = new CGnuplot();

	// On masque la communication
	gpl->setLogging(false);
	
	// Le signal est défini apres le fork, pour qu'il ne concerne que ce
	// process, et pas gnuplot.
	signal(SIGINT,sighdl);

	// Definit des options de gnuplot et une fonction
	gpl->plot("set terminal png;");
	gpl->plot("set surface;set hidden3d;set isosamples 15,15");
	gpl->plot("set mouse");
	gpl->plot("f(a,b)=((a*a)>(b*b))?a:b");


	for (i=1;i<100;i+=10)
	{
		// Affichage de la fonction. La fonction rend la main quand gnuplot
		// rend la main...
		gpl->plot("set output \"film%04d.png\"",i);
		gpl->plot("splot [-4.5:4.5][-4.5:4.5][-1:1] exp(-%f*abs(f(x,y)))*cos(%f*f(x,y))",0.05*abs(50-(i%100)),1+0.005*abs(50-(i%100)));

		if (signal_received>=1) break;
	}

	// On suppose que gnuplot sait quitter rapidement et proprement
	gpl->setDieTime(0);
	delete gpl;
	
	return 0;
}
