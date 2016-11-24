/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#ifndef CGNUPLOT_H
#define CGNUPLOT_H

#include "CShell.h"

namespace cgnuplot {
    /**
     * \class CGnuplot : 
     * permet de controler gnuplot
     * **/
    class CGnuplot : public CShell
    {
        protected :
            unsigned int dietime_sec;
            float timeout;

        public :
            /** 
             * Constructeur : instancie gnuplot et attend le prompt.
             * Prerequis : gnuplot doit etre dans le PATH !
             * binname : le nom (avec eventuellement le path) de gnuplot
             * options : options a passer a gnuplot.
             * */
            CGnuplot(const char * binname = "/usr/bin/gnuplot", const char * options = "-noraise");

            /**
             * Veut ont afficher la communication avec gnuplot ? 
             * */
            void setLogging(bool l);
            /**
             * Combien de temps faut il laisser a gnuplot pour sortir apres exit.
             * Par defaut : 0 sec (3.8j). Il faut 1 sec pour la 3.7
             * */
            void setDieTime(unsigned int t_sec);

            void setTimeout(float tout) {timeout = tout;}
            float getTimeout() const {return timeout;}
            /**
             * Destructeur : termine gnuplot
             * */
            ~CGnuplot();

            /**
             * Interromp l'execution de la commande en cours 
             * */
            void interrupt() {sendsig(SIGINT);}


            /**
             * Envoie la commande cmd a gnuplot et rend la main quand gnuplot
             * affiche son prompt.
             * Meme syntaxe que printf. Par exemple : 
             * send("plot [-3.14:3.14] %s(x/ %d) w l","sin", 4.0);
             * **/
            void plot(const char * cmd, ...);

    };
};

#endif // CGNUPLOT_H
