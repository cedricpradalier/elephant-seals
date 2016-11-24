/**Signature>
* Author      : Cedric Pradalier 
* Universite  : INRIA - GRAVIR - INPG
* Email       : cedric.pradalier@inrialpes.fr
* Contexte    : These MESR 
* Date        : 2001 - 2004
* License     : Libre (???)
<Signature**/
#include <sched.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>
#include <signal.h>
#include <sys/wait.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <termios.h>
#include "cgnuplot/CShell.h"

// #define DEBUG
using namespace cgnuplot;

static int waitdata(int fd,unsigned int tsec, unsigned int tmsec)
{
	fd_set rfs;
	struct timeval to;
	to.tv_sec = tsec;
	to.tv_usec = tmsec*1000;
	FD_ZERO(&rfs);FD_SET(fd,&rfs);
#ifdef DEBUG
	// printf("Waiting data\n");fflush(stdout);
#endif
	return select(fd+1,&rfs,NULL,NULL,&to);
}

CShell::CShell()
{
	echo = log = false;
	p_in[0] = p_in[1] = -1;
	p_out[0] = p_out[1] = -1;
	pid = 0;
}

void CShell::sendsig(int signum)
{
	kill(pid,signum);
}

bool CShell::start(const char * command) 
{
	int fl,r;
	
	r = pipe(p_in);
    if (r < 0) {
        perror("Creating input pipe: ");
        return false;
    }
	// printf("Pipe in: r %d [%d %d]\n",r,p_in[0],p_in[1]);
	r = pipe(p_out);
    if (r < 0) {
        perror("Creating output pipe: ");
        return false;
    }
	// printf("Pipe out: r %d [%d %d]\n",r,p_out[0],p_out[1]);
	// fl = fcntl(p_in[0],F_GETFL); fcntl(p_in[0],F_SETFL, fl | O_SYNC);
	// fl = fcntl(p_out[1],F_GETFL); fcntl(p_out[1],F_SETFL, fl | O_SYNC);
	pid = fork();
	if (pid < 0) {
		perror("Forking new process: ");
	}
	if (pid == 0)
	{
		// fils
#ifdef DEBUG
		printf("Fils : executing : \"sh -c %s\"\n", command);
#endif
		setenv("TERM","vt100",1);
		close(p_out[0]);
		close(p_in[1]);
		dup2(p_in[0],STDIN_FILENO);
		dup2(p_out[1],STDOUT_FILENO);
#if 1
		if (execlp("sh","sh","-c",command,NULL)) {
			perror("Execlp failed");
			return false;
		}
#else
		perror("Testing");
		if (execlp("ping","ping","localhost",NULL)) {
			perror("Execlp failed");
			return false;
		}
#endif
		// never return
	}
	// pere
	close(p_out[1]);
	close(p_in[0]);

	fl = fcntl(p_in[1],F_GETFL); fcntl(p_in[1],F_SETFL, fl | O_SYNC);
	fl = fcntl(p_out[0],F_GETFL); fcntl(p_out[0],F_SETFL, fl | O_NONBLOCK );

	fpin = fdopen(p_in[1],"w");

	return true;
}

bool CShell::terminate(bool killhard)
{
	if (pid == 0) return false;
	flushOutput(1.0);
	close(p_out[0]); p_out[0] = 0;
	close(p_in[1]); p_in[1] = 0;

	usleep(10000);
	
	kill(pid,SIGTERM);
	if (killhard)
	{
		//printf("Waiting to Kill\n");
		sleep(1);
		kill(pid,SIGKILL);
	}
	pid = 0;
	
	return true;
}


CShell::~CShell()
{
	terminate(true);
}

void CShell::send(const char * cmd, ...)
{
	va_list args;
	va_start(args, cmd);
	if (echo) {vprintf(cmd,args);fflush(stdout);}
	va_end(args);
	va_start(args, cmd);
	vfprintf(fpin,cmd,args);
	va_end(args);
	fflush(fpin);
	sched_yield();
}

void CShell::flushOutput(double timeout_sec)
{
	int n;
	char c[8192];
	double t0,t1;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	t0 = tv.tv_sec+(1e-6*tv.tv_usec);
	
	while (1)
	{
		n = read(p_out[0],c,8192);
		if (n <= 0) return ;
		if (log) {fwrite(c,1,n,stdout);}

		gettimeofday(&tv,NULL);
		t1 = tv.tv_sec+(1e-6*tv.tv_usec);
		if (t1-t0>timeout_sec) return;
	}
}

bool CShell::waitPrompt(const char * prompt, 
		double timeout_sec)
{
	unsigned int n = strlen(prompt);
	char buffer[n];
	unsigned int i,w=0;
	double t0,t1;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	t0 = tv.tv_sec+(1e-6*tv.tv_usec);

	while (1)
	{
		int r =  waitdata(p_out[0],0,10);
		gettimeofday(&tv,NULL);
		t1 = tv.tv_sec+(1e-6*tv.tv_usec);
		if (t1-t0>timeout_sec) return false;
		if (r <= 0) continue;

		// WaitData dit qu'il y a des donnees
		assert (read(p_out[0],&(buffer[w]),1) == 1);
		if (log) { fprintf(stdout,"%c",buffer[w]); }

		// Il faudrait un automate !!!
		w += 1;
		if (w>=n)
		{
			if (strncmp(buffer,prompt,n)==0)
				return true;
			for (i=0;i<n-1;i++)
				buffer[i]=buffer[i+1];
			w -= 1;
		}
	}
}

