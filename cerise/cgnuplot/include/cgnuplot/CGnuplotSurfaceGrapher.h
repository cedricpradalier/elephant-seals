/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
/// Ronan Le Hy, 2004

#ifndef CGNUPLOT_SURFACE_GRAPHER_H
#define CGNUPLOT_SURFACE_GRAPHER_H

#include "CGnuplot.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <string>
using namespace std;

namespace cgnuplot {
    /**
     * Manages drawing 3D surfaces in a CGnuplot object.
     **/
    class CGnuplotSurfaceGrapher
    {
        public:
            CGnuplotSurfaceGrapher(CGnuplot& cgnuplot);
            virtual ~CGnuplotSurfaceGrapher();

            void init_gnuplot(const int& x_samples = 10, 
                    const int& y_samples = 10, 
                    const bool& mouse = true);

            /**
             * Draws a z=f(x,y) surface w here v[i*width+j] contains f(i,j).
             * @param v the vector containing data values
             * @param width width of the data matrix stored in the vector
             * @param hheight height of the data matrix stored in the matrix, default is -1: calculate from
             *        vector size and width
             * @param windowNumber an integer identifying the x11 window in which the plot will be done 
             *        (default is -1: let gnuplot use the last one)
             * @param title the title of the plot -- must not contain double quotes ("), or else, or else what? exactly!
             * @return true on success, false on failure (caution: very lightly tested)
             **/
            template<typename T>
                bool
                do_draw(const vector<T>& v, 
                        const int& width, 
                        const int& hheight = -1, 
                        const string& title = "")
                {
                    if (width == 0 || v.size() == 0)
                    {
                        /* 		cerr << "draw got width=0 -- not drawing" << endl; */
                        return false;
                    }

                    int height = hheight;
                    if (height == -1)
                    {
                        height = v.size() / width;
                    }

                    if (title != "")
                    {
                        _cgnuplot.send("splot '-' title \"%s\" with lines\n", title.c_str());		
                    }
                    else
                    {
                        _cgnuplot.send("splot '-' with lines\n");
                    }

                    int i = 0;
                    for (int w=0; w<width; ++w)
                    {
                        for (int h=0; h<height; ++h)
                        {
                            _cgnuplot.flushOutput(0.01);
                            _cgnuplot.send("%d %d %g\n", w, h, double(v[i]));
                            ++i;
                        }
                        _cgnuplot.send("\n");
                    }
                    _cgnuplot.send("e\n");
                    _cgnuplot.send("print \"ready\"\n");

                    return _cgnuplot.waitPrompt("ready", _cgnuplot.getTimeout());
                };

            template<typename T>
                bool
                draw(const vector<T>& v, 
                        const int& width, 
                        const int& hheight = -1, 
                        int windowNumber = -1, 
                        const string& title = "")
                {
                    if (windowNumber != -1)
                    {
                        _cgnuplot.send("set terminal x11 %d\n", windowNumber);
                    }

                    return do_draw(v, width, hheight, title);
                };


            template<typename T>
                bool
                draw_to_file(const string& filename,
                        const vector<T>& v, 
                        const int& width, 
                        const int& hheight = -1, 
                        const string& title = "")
                {
                    _cgnuplot.send("set terminal x11\n");

                    _cgnuplot.send("set terminal png medium\n");
                    _cgnuplot.send("set output \"%s\"\n", filename.c_str());

                    return do_draw(v, width, hheight, title);
                };

        protected:
            CGnuplot& _cgnuplot;
    };
};

#endif /* CGNUPLOT_SURFACE_GRAPHER_H */
