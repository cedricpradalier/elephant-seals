/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#include "cgnuplot/CGnuplot.h"
#include "cgnuplot/CGnuplotSurfaceGrapher.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ext/stdio_filebuf.h>

using namespace std;
using namespace cgnuplot;

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

CGnuplotSurfaceGrapher::CGnuplotSurfaceGrapher(CGnuplot& cgnuplot) :
    _cgnuplot(cgnuplot)
{
    // empty
}

CGnuplotSurfaceGrapher::~CGnuplotSurfaceGrapher()
{
}

/**
 * Non mandatory. You can also initialise the plot as you want with send() and plot() on the CGnuplot object.
 **/
void CGnuplotSurfaceGrapher::init_gnuplot(const int& x_samples, 
					  const int& y_samples, 
					  const bool& mouse)
{
    _cgnuplot.plot("set surface; set hidden3d;set isosamples %d,%d", x_samples, y_samples);

    if (mouse)
    {
	_cgnuplot.plot("set mouse");
    }
}
