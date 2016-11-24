/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#include <sched.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include "cgnuplot/CGnuplot.h"

using namespace cgnuplot;

CGnuplot::CGnuplot(const char * binname, const char * options) : CShell()
{
	char command[8192];
	dietime_sec = 1;
	sprintf(command,"%s %s 2>&1",binname,options);
	assert(start(command));
	setTimeout(3.0);
	setLog(false);
	setEcho(false);
	send("set terminal x11;print \"ready\"\n");
	assert(waitPrompt("ready",3));
}

CGnuplot::~CGnuplot()
{
	CShell::send("exit\n");
	if (dietime_sec>0)
	{
		//old gnuplot version
		sleep(dietime_sec);
	}
	terminate(false);
}

void CGnuplot::setLogging(bool l)
{
	CShell::setLog(l);
	CShell::setEcho(l);
}

void CGnuplot::setDieTime(unsigned int t_sec)
{
	dietime_sec = t_sec;
}

void CGnuplot::plot(const char * command, ...)
{
	char cmd[8192];
	va_list args;
	va_start(args, command);
	vsprintf(cmd,command,args);
	va_end(args);
	strcat(cmd,"\n");

	send(cmd);
	send("print \"ready\"\n");

	waitPrompt("ready", timeout);
}


