/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#ifndef CSHELL_H
#define CSHELL_H

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

namespace cgnuplot {
    /**
     * \class CShell : 
     * permet de controler un programme se comportant comme un shell,
     * TODO : ameliorer la gestion d'erreur...
     * **/
    class CShell 
    {
        protected : 
            int p_in[2];
            int p_out[2];
            FILE *fpin;
            pid_t pid;
            bool log;
            bool echo;

        public :
            /** 
             * Constructeur : ne fait rien !
             * */
            CShell();

            /**
             * Veut-on visualiser les sorties, pendant flushOutput ou waitPrompt ?
             * si l == true, on les affiche
             * **/
            void setLog(bool l) {log = l;}
            /**
             * Veut-on afficher les commandes envoyées par send ?
             * si e == true, on les affiche
             * **/
            void setEcho(bool e) {echo = e;}

            /**
             * execute le programme command dans sh: fork + pipe + execlp (cf. man)
             * par exemple :
             * start("gnuplot --noraise -");
             * */
            bool start(const char * command);

            /**
             * Termine le programme commandé (SIGTERM) et ferme les pipe
             * Si killhard est vrai, laisse une seconde pour terminer et on envoie
             * un SIGKILL !
             * */
            bool terminate(bool killhard=true);

            /**
             * Envoie un signal au programme commandé.
             * Attention, si celui-ci ne sait pas le recevoir,
             * il risque de mourrir et le prochain send/flush
             * provoquera un sigpipe
             * */
            void sendsig(int signum=SIGINT);

            /**
             * Destructeur : termine le programme controllé
             * */
            ~CShell();

            /**
             * Envoie la commande cmd au programme controllé
             * Meme syntaxe que printf. Par exemple : 
             * send("plot [-3.14:3.14] %s(x/ %d) w l","sin", 4.0);
             * **/
            void send(const char * cmd, ...);

            /**
             * Consomme les données sur la sortie standard du programme controllé,
             * pendant au plus timeout_sec secondes.
             * */
            void flushOutput(double timeout_sec);

            /**
             * Scan la  sortie standard du programme controllé,
             * a la recherche de prompt (TODO : utiliser un automate)
             * pendant au plus timeout_sec seconde.
             * Renvoie false si la sortie est due au timeout.
             * */
            bool waitPrompt(const char * prompt, double timeout_sec);
    };
};
#endif // CSHELL_H
