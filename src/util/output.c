/***************************************************************
FronTier is a set of libraries that implements differnt types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

/*
*		Routines for Random Access Output:
*
*	The routine foutput() allows programs to maintain output files
*	as linked lists of variable length blocks.   To start a new
*	block, just call foutput().   foutput() returns 1, or 0 on error
*	which generally indicates that the output was not to a
*	seekable device.  output() is a special case of foutput() where
*	the output file is simply standard output.
*
*	The corresponding routine next_output() skips to the beginning
*	of the next output block in a file that was previously written with
*	foutput() and returns 1, or 0 on error.   On other files nothing
*	is done.  Similarly prev_output() skips to the previous output
*	block, assuming next_output() was already called.
*
*	A file is recognized as having been generated by
*	foutput() if it starts with newline in bytes 1, 22,  # in bytes 2,23
*	and an integer in bytes 3-11, space in byte 12, and another integer
*	in bytes 13-21.   The routine check_output() may be called to perform
*	this test.   The routines is_start_output() and is_end_output()
*	return 1 if there is no previous or later output mark in the file,
*	or 0 otherwise.  They do not actually read or move the file pointer.
*    
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cdecs.h>
#include <vmalloc.h>

#undef	fclose

struct _IOUTPUT { 
	IO_TYPE Io_type;
	long    save;
	boolean    output_blocked;
	boolean    check_for_blocked_output;
	long    next, current, prev;
};
typedef struct _IOUTPUT IOUTPUT;

struct _READ_LIST {
	IOUTPUT Output;
	struct _READ_LIST *prev, *next;
};
typedef struct _READ_LIST READ_LIST;

struct _WRITE_LIST {
	FILE *file;
	boolean seekable;
	long previous, last;
	struct _WRITE_LIST *prev, *next;
};
typedef struct _WRITE_LIST WRITE_LIST;


enum _ADD_OR_DELETE {
	DELETE = 0,
	ADD    = 1
};
typedef enum _ADD_OR_DELETE ADD_OR_DELETE;

	/* LOCAL function Prototypes */
LOCAL	IOUTPUT	*index_of_read_file(FILE*,ADD_OR_DELETE);
LOCAL	WRITE_LIST *index_of_write_file(FILE*,ADD_OR_DELETE);
LOCAL	boolean	reopen_foutput(WRITE_LIST*);
LOCAL	void	set_new_read_file(FILE*,IOUTPUT*);

LOCAL	boolean	determine_io_type_blocked = NO;
LOCAL	size_t    default_read_float_size = sizeof(double);
LOCAL	boolean      default_reverse_endian = NO;
LOCAL	FT_ENDIAN default_read_endian = FT_UNKNOWN_ENDIAN;


EXPORT boolean output(void)
{
	return foutput(stdout);
}		/*end output*/


EXPORT boolean foutput(
	FILE		*file)
{
	WRITE_LIST	*wlist;
	long		current;

	wlist = index_of_write_file(file,ADD);

	if (!wlist->seekable)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in foutput(), "
	                    "file not seekable\n");
	    return NO;
	}

	current = ftell(file);

	if (fseek(file,wlist->last,SEEK_SET) < 0)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in append_output(), "
	                    "fseek(file,0L,SEEK_END) < 0\n");
	    return NO;
	}

	(void) fprintf(file,"\n#%9ld",current);
	(void) fseek(file,current,SEEK_SET);
	(void) fprintf(file,"\n#%9ld %9ld\n#",-1L,wlist->previous);

	wlist->last = current;
	wlist->previous = current;
	return YES;
}		/*end foutput*/


EXPORT	boolean append_output(
	FILE	*file)
{
	WRITE_LIST	*wlist;
	long		current;

	if (!check_output(file))
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in append_output(), "
	                    "check_output() returns NO\n");
	    return NO;
	}

	/*Position at end of output blocks*/

	(void) rewind_read_file(file,NULL);
	while (next_output(file));
	wlist = index_of_write_file(file,ADD);
	(void) fseek(file,-23L,SEEK_CUR);
	current = ftell(file);
	wlist->last = current;
	wlist->previous = current;
	if (fseek(file,0L,SEEK_END) < 0)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in append_output(), "
	                    "fseek(file,0L,SEEK_END) < 0\n");
	    return NO;
	}
	return YES;
}		/*end append_output*/


EXPORT	boolean	erase_last_foutput(
	FILE	*file)
{
	WRITE_LIST	*wlist;
	long		current;
	long		return_to;

	return_to = ftell(file);
	if (!check_output(file))
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in erase_last_foutput(), "
	                    "check_output() returns NO\n");
	    return NO;
	}
	(void) rewind_read_file(file,NULL);
	while (next_output(file));
	(void) prev_output(file);
	wlist = index_of_write_file(file,ADD);
	(void)fseek(file,-23L,SEEK_CUR);
	current = ftell(file);
	wlist->last = current;
	wlist->previous = current;
	(void) fprintf(file,"\n#%9ld",-1L);
	(void)fseek(file,return_to,SEEK_SET);
	return YES;
}		/*end erase_last_foutput*/


EXPORT int Fclose(
	FILE		*file)
{
	(void) index_of_read_file(file,DELETE);
	(void) index_of_write_file(file,DELETE);
	return fclose(file);
}		/*end Fclose*/


LOCAL WRITE_LIST *index_of_write_file(
	FILE		*file,
	ADD_OR_DELETE	flag)
{
	WRITE_LIST	*wlist;
	struct stat	x;
	static WRITE_LIST *head_of_write_list = NULL,
	                  *tail_of_write_list = NULL;

	for (wlist = head_of_write_list; wlist != NULL; wlist = wlist->next)
	    if (file == wlist->file) break;

	if (flag == DELETE)
	{
	    if (wlist == NULL)
	        return NULL;

	    if (wlist == head_of_write_list)
	        head_of_write_list = wlist->next;
	    if (wlist == tail_of_write_list)
	        tail_of_write_list = wlist->prev;
	    if (wlist->prev)
	        wlist->prev->next = wlist->next;
	    if (wlist->next)
	        wlist->next->prev = wlist->prev;
	    free(wlist);
	    return NULL;
	}

	if (wlist != NULL)
	    return wlist;

	scalar(&wlist,sizeof(WRITE_LIST));
	if (wlist == NULL)
	    return NULL;
	wlist->next = NULL;
	wlist->prev = tail_of_write_list;
	if (head_of_write_list == NULL)
	    head_of_write_list = wlist;
	else
	    tail_of_write_list->next = wlist;
	tail_of_write_list = wlist;

	wlist->file = file;
	wlist->previous = -1L;
	wlist->last =  0L;
	wlist->seekable = ( (fstat(fileno(file),&x)==0) &&
	                        (x.st_nlink != 0 /*!pipe*/) && 
	                        !isatty(fileno(file)) ) ? YES : NO;
	wlist->seekable = reopen_foutput(wlist);
	return wlist;
}		/*end index_of_write_file*/


LOCAL	boolean reopen_foutput(
	WRITE_LIST	*wlist)
{
	FILE		*file;
	long		save;
	int		flag;

	if (!wlist->seekable)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in reopen_foutput, "
	                    "wlist->seekable = NO\n");
	    return NO;
	}

	file = wlist->file;
	save = ftell(file);
	if (save == 0)
	    return YES;
	flag = fcntl(fileno(file),F_GETFL,0);
	if ((flag&O_RDWR) == 0      	  ||
	    (!check_output(file))	  ||
	    (fseek(file,0L,SEEK_END) < 0) ||
	    (index_of_read_file(file,ADD) == NULL))
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in reopen_foutput, "
	                    "searchable tests fail\n");
	    return NO;
	}
	while (next_output(file));
	wlist->previous = wlist->last = ftell(file)-23;
	if (fseek(file,save,SEEK_SET) < 0)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in reopen_foutput(), "
	                    "fseek(file,0L,SEEK_END) < 0\n");
	    return NO;
	}
	return YES;
}		/*end reopen_foutput*/


LOCAL IOUTPUT *index_of_read_file(
	FILE		*file,
	ADD_OR_DELETE	flag)
{
	READ_LIST	*rlist;
	static READ_LIST *head_of_read_list = NULL,
	                 *tail_of_read_list = NULL;

	debug_print("foutput","Entered index_of_read_file()\n");
	if (debugging("foutput"))
	{
	    (void) printf("searching for file %p\n",file);
	    (void) printf("flag = %s\n",
	                  (flag == DELETE) ? "DELETE" : "ADD");
	}

	for (rlist = head_of_read_list; rlist != NULL; rlist = rlist->next)
	    if (file == rlist->Output.Io_type.file)
	        break;

	if (debugging("foutput"))
	{
	    if (rlist != NULL)
	        (void) printf("Entry IOUTPUT %p found\n",&rlist->Output);
	}

	if (flag == DELETE)
	{
	    if (rlist != NULL)
	    {
	        if (rlist == head_of_read_list)
	            head_of_read_list = rlist->next;
	        if (rlist == tail_of_read_list)
	            tail_of_read_list = rlist->prev;
	        if (rlist->prev)
	            rlist->prev->next = rlist->next;
	        if (rlist->next)
	            rlist->next->prev = rlist->prev;
	        free(rlist);
	    }
	    debug_print("foutput","Left index_of_read_file()\n");
	    return NULL;
	}

	if (rlist != NULL)
	{
	    if (debugging("foutput"))
	        (void) printf("check_for_blocked_output = %s\n",
		              y_or_n(rlist->Output.check_for_blocked_output));
	    debug_print("foutput","Left index_of_read_file()\n");
	    return &rlist->Output;
	}

	/* A new file */

	if (debugging("foutput"))
	    (void) printf("Adding new entry\n");

	scalar(&rlist,sizeof(READ_LIST));
	if (rlist == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in index_of_read_file(), "
	                    "can't allocate rlist\n");
	    debug_print("foutput","Left index_of_read_file()\n");
	    return NULL;
	}
	rlist->next = NULL;
	rlist->prev = tail_of_read_list;
	if (head_of_read_list == NULL)
	    head_of_read_list = rlist;
	else
	    tail_of_read_list->next = rlist;
	tail_of_read_list = rlist;
	set_new_read_file(file,&rlist->Output);
	if (debugging("foutput"))
	    (void) printf("Entry IOUTPUT %p added\n",&rlist->Output);
	debug_print("foutput","Left index_of_read_file()\n");
	return &rlist->Output;
}		/*end index_of_read_file*/


LOCAL	void	set_new_read_file(
	FILE		*file,
	IOUTPUT		*oput)
{
	oput->Io_type.file = file;
	oput->Io_type.read_float_size = default_read_float_size;
	oput->Io_type.read_endian     = default_read_endian;
	oput->Io_type.reverse_endian  = default_reverse_endian;
	oput->Io_type.cpu_float_size  = sizeof(double);
	oput->Io_type.cpu_endian = ft_endian_type();
	oput->next = 0L;
	oput->prev = -1L;
	oput->current = 0L;
	oput->save = 0L;
	oput->output_blocked = NO;
	oput->check_for_blocked_output = YES;
}		/*end set_new_read_file*/


EXPORT boolean next_output(
	FILE		*file)
{
	IOUTPUT		*oput;
	char		s[24];
	boolean		KeepReading = YES;

	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in next_output(), "
	                    "index_of_read_file() returns NULL\n");
	    return NO;
	}
	oput->current = ftell(file);
	while (KeepReading == YES)
	{
	    if (oput->next < 0) /*End of file reached*/
	        return NO;
	    if ((fseek(file,oput->next,SEEK_SET) == -1) ||
		(fread((void*)s,1,23,file) != 23)       ||
		(s[0]!='\n') 			        || 
		(s[21]!='\n') 			        || 
		(s[1]!='#') 			        || 
		(s[22]!='#'))
	    {
	        (void) printf("WARNING in next_output(), "
	                      "can't find next output block, "
			      "oput->next = %ld\n",oput->next);
	        return NO;
	    }
	    oput->next = atol(s+2);
	    oput->prev = atol(s+12);
	    if (ftell(file) > oput->current)
	        break;
	}
	return YES;
}		/*end next_output*/


EXPORT boolean is_end_output(
	FILE		*file)
{
	IOUTPUT		*oput;

	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in is_end_output(), "
	                    "index_of_read_file() returns NULL\n");
	    return NO;
	}
	return (oput->next < 0) ? YES : NO;
}		/*end is_end_output*/


EXPORT boolean is_start_output(
	FILE		*file)
{
	IOUTPUT		*oput;

	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in is_start_output(), "
	                    "index_of_read_file() returns NULL\n");
	    return NO;
	}
	return (oput->prev < 0) ? YES : NO;
}		/*end is_start_output*/


EXPORT boolean prev_output(
	FILE		*file)
{
	IOUTPUT		*oput;
	char		s[24];
	boolean		KeepReading = YES;

	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in prev_output(), "
	                    "index_of_read_file() returns NULL\n");
	    return NO;
	}
	oput->current = ftell(file);
	while (KeepReading == YES)
	{
	    if (oput->prev < 0) /* Beginning of file reached */
	        return NO;
	    if ((fseek(file,oput->prev,SEEK_SET) == -1) ||
		(fread((void*)s,1,23,file)!=23)		||
		(s[0]!='\n')	 			|| 
		(s[21]!='\n') 				|| 
		(s[1]!='#') 				|| 
		(s[22]!='#'))
	    {
	        (void) printf("WARNING in prev_output(), "
	                      "can't find prev output block "
			      "oput->prev = %ld\n",oput->prev);
	        return NO;
	    }
	    oput->next = atol(s+2);
	    oput->prev = atol(s+12);
	    if (ftell(file) < oput->current)
	        break;
	}
	return YES;
}		/*end prev_output*/


EXPORT boolean check_output(
	FILE		*file)
{
	IOUTPUT		*oput;
	int		n;
	char		s[24];
	static long	save,next,prev;
	IMPORT boolean	fgetstring_debug;

	debug_print("foutput","Entered check_output(), file = %p\n",file);
	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in check_output(), "
	                    "index_of_read_file() returns NULL\n");
	    debug_print("foutput","Left check_output(), ans = %s\n",y_or_n(NO));
	    return NO;
	}
	if (oput->check_for_blocked_output)
	{
	    if (debugging("foutput"))
	        (void) printf("Testing file %p for blocked output\n",file);
	    oput->output_blocked = NO;
	    save = ftell(file);
	    if (fseek(file,0L,SEEK_SET) < 0)
	    {
	        if (debugging("foutput"))
	          (void) printf("WARNING in check_output(), "
	                        "fseek(file,0L,SEEK_SET) < 0\n");
	        debug_print("foutput","Left check_output(), ans = %s\n",y_or_n(NO));
	        return NO;
	    }
	    for (n=0; n<23; ++n)
	        s[n] = getc(file);
	    s[23] = '\0';
	    oput->output_blocked = YES;
	    if (s[0]!='\n' || s[21]!='\n' || s[1]!='#' || s[22]!='#')
	    {
	        if (debugging("foutput"))
	          (void) printf("WARNING in check_output(), "
	                        "Can't find leading block printout\n");
	        oput->output_blocked = NO;
	    }
	    if (sscanf(s+2,"%ld %ld",&next,&prev)!=2)
	    {
	        if (debugging("foutput"))
	          (void) printf("WARNING in check_output(), "
	                        "Can't read prev and next addresses\n");
	        oput->output_blocked = NO;
	    }

	    if (fgetstring_debug || debugging("foutput"))
	    {
	        if (oput->output_blocked)
	            (void) fprintf(stderr,"Using output trace in file\n");
	        else
	            (void) fprintf(stderr,"No output trace in file\n");
	    }
	    oput->check_for_blocked_output = NO;
	    if (fseek(file,save,SEEK_SET)<0)
	    {
	        if (debugging("foutput"))
	          (void) printf("WARNING in check_output(), "
	                        "fseek(file,save,SEEK_SET) < 0\n");
	        oput->output_blocked = NO;
	        debug_print("foutput","Left check_output(), ans = %s\n",y_or_n(NO));
	        return NO;
	    }
	    if (!determine_io_type_blocked)
	    {
	        if (debugging("foutput"))
	            (void) printf("Calling determine_io_type()\n");
	        determine_io_type(file,&oput->Io_type);
	    }
	}
	if (debugging("foutput") && !oput->output_blocked)
	{
	    (void) printf("WARNING in check_output(), "
	                  "!oput->output_blocked\n");
	}
	debug_print("foutput","Left check_output(), ans = %s\n",
	      y_or_n(oput->output_blocked));
	return oput->output_blocked;
}		/*end check_output*/


EXPORT void trace_foutput(
	FILE	*file)
{
	static const char *TRACE_MARK = "#$$$002$$$";	/* Version Number */
	if (foutput(file))
	    (void) fprintf(file,"\n%s\n",TRACE_MARK);
}		/*end trace_foutput*/


EXPORT void trace_output(void)
{
	trace_foutput(stdout);
}		/*end trace_output*/


EXPORT const char *next_output_line_containing_string(
	FILE	   *file,
	const char *strng)
{
	IMPORT       boolean fgetstring_debug;
	static const int  MAX_LINE_LEN = 2048;
	IOUTPUT	          *oput;
	boolean	          blank;
	static char       *line = NULL;

	if (fgetstring_debug)
	    (void) fprintf(stderr,"next_output_line_containing_string(%s) ",strng);
	if (line == NULL)
	    uni_array(&line,MAX_LINE_LEN+2,sizeof(char));

	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (fgetstring_debug)
	        (void) fprintf(stderr,"FAIL\n");
	    return NULL;
	}

	if (oput->check_for_blocked_output)
	{
	    if (fgetstring_debug)
	        (void) fprintf(stderr,"call check_output() ");
	    (void) check_output(file); 
	}

	switch (strng[0])
	{
	case '\0':
	case '\n':
	    blank = YES;
	    break;
	default:
	    blank = NO;
	    break;
	}

	if (oput->output_blocked)
	{
	    static OUTPUT *oput = NULL;
	    oput = save_read_file_variables(file,oput);
	    while (next_output(file))
	    {
	        if (blank == YES)
		{
	            (void) fgets(line,MAX_LINE_LEN,file);
		}
	        else
	        {
	            boolean KeepReading = YES;

	            while (KeepReading == YES)
	            {
	                if (fgets(line,MAX_LINE_LEN,file) == NULL)
	                {
	                    reset_read_file_variables(oput);
	                    if (fgetstring_debug)
	                        (void) fprintf(stderr,"FAIL\n");
	                    return NULL;
	                }
	                if (line[0]!='\n')
	                    break;
	            }
	        }
	        if (sgetstring(line,strng) != NULL)
	            return line;
	    }
	    reset_read_file_variables(oput);
	    if (fgetstring_debug)
	        (void) fprintf(stderr,"FAIL\n");
	    return NULL;
	}
	else
	{
	    long save = ftell(file);
	    while (fgets(line,MAX_LINE_LEN,file) != NULL) 
	    {
	        if (sgetstring(line,strng) != NULL)
	            return line;
	    }
	    (void) fseek(file,save,SEEK_SET);
	    if (fgetstring_debug)
	        (void) fprintf(stderr,"FAIL\n");
	    return NULL;
	}
}		/*end next_output_line_containing_string*/


EXPORT OUTPUT *save_read_file_variables(
	FILE		*file,
	OUTPUT		*savoput)
{
	IOUTPUT		*oput;
	IOUTPUT		*soput;

	debug_print("foutput","Entered save_read_file_variables()\n");
	if (file == NULL)
	{
	    debug_print("foutput","Left save_read_file_variables()\n");
	    return savoput;
	}
	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in save_read_file_variables(), "
	                    "index_of_read_file() returns NULL\n");
	    debug_print("foutput","Left save_read_file_variables()\n");
	    return savoput;
	}
	soput = (IOUTPUT*)savoput;
	if (soput == NULL)
	  scalar(&soput,sizeof(IOUTPUT));
	*soput = *oput;
	soput->save = ftell(file);
	debug_print("foutput","Left save_read_file_variables()\n");
	return (OUTPUT*)soput;
}		/*end save_read_file_variables*/


EXPORT boolean rewind_read_file(
	FILE		*file,
	OUTPUT		*rwnd)
{
	IOUTPUT		*oput;

	debug_print("foutput","Entered rewind_read_file()\n");
	if ((oput = index_of_read_file(file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in rewind_read_file(), "
	                    "index_of_read_file() returns NULL\n");
	    debug_print("foutput","Left rewind_read_file()\n");
	    return NO;
	}
	if (rwnd != NULL)
	    (void) save_read_file_variables(file,rwnd);
	oput->next = 0L;
	oput->current = 0L;
	oput->prev = -1L;
	rewind(file);
	oput->save = ftell(file);
	debug_print("foutput","Left rewind_read_file()\n");
	return YES;
}		/*end rewind_read_file*/


EXPORT void reset_read_file_variables(
	OUTPUT		*rset)
{
	IOUTPUT	*oput;
	IOUTPUT *reset_oput;
	boolean    check_for_blocked_output;
	boolean    output_blocked;

	debug_print("foutput","Entered reset_read_file_variables()\n");
	reset_oput = (IOUTPUT*)rset;
	if ((oput = index_of_read_file(reset_oput->Io_type.file,ADD)) == NULL)
	{
	    if (debugging("foutput"))
	      (void) printf("WARNING in reset_read_file_variables(), "
	                    "index_of_read_file() returns NULL\n");
	    debug_print("foutput","Left reset_read_file_variables()\n");
	    return;
	}
	(void) fseek(reset_oput->Io_type.file,reset_oput->save,SEEK_SET);
	output_blocked = NO;
	if ((oput->check_for_blocked_output == NO) ||
	    (reset_oput->check_for_blocked_output == NO))
	{
	    if (oput->output_blocked || reset_oput->output_blocked)
	        output_blocked = YES;
	    check_for_blocked_output = NO;
	}
	else
	    check_for_blocked_output = YES;
	*oput = *reset_oput;
	oput->check_for_blocked_output = check_for_blocked_output;
	oput->output_blocked = output_blocked;
	debug_print("foutput","Left reset_read_file_variables()\n");
	return;
}		/*end reset_read_file_variables*/


EXPORT	const IO_TYPE *open_close_file(
	const char * const *nfiles,
	int	   i_open,
	int	   i_close,
	int	   reopen_stdin)
{
	static int	   nf;
	static IOUTPUT	   *olist = NULL;
	static const char * const *currlist = NULL;

	debug_print("foutput","Entered open_close_file()\n");
	if ((currlist != nfiles) || (nfiles == NULL))
	{
	    int		i;

	    /* Destroy any prexisting list */
	    if (olist != NULL)
	    {
	        for (i = 0; i < nf; ++i)
	        {
	            if (olist[i].Io_type.file != NULL)
	                (void) Fclose(olist[i].Io_type.file);
	        }
	        free(olist);
	        olist = NULL;
	    }
	    if (nfiles == NULL)
	    {
	        debug_print("foutput","Left open_close_file()\n");
	        return NULL;
	    }

	    /* Count current files and set up new file control structures*/
	    for (nf = 0, currlist = nfiles; *nfiles; ++nf, ++nfiles);
	    if (nf > 0)			
	        uni_array(&olist,nf,sizeof(IOUTPUT));
	    for (i = 0; i < nf; ++i)
	        set_new_read_file(NULL,olist+i);
	}

	if (*currlist == NULL)
	{
	    static IO_TYPE Io_type;
	    determine_io_type(stdin,&Io_type);
	    debug_print("foutput","Left open_close_file()\n");
	    return &Io_type;
	}

	/* closing part */
	if (i_close >= 0 && i_close<nf && olist[i_close].Io_type.file!=NULL)
	{
	    (void) save_read_file_variables(olist[i_close].Io_type.file,
	                                    &olist[i_close]);
	    (void) fclose(olist[i_close].Io_type.file);
	    olist[i_close].Io_type.file = NULL;
	}

	/* opening part */
	if (i_open >= 0 && i_open < nf)
	{
	    olist[i_open].Io_type.file = (reopen_stdin == YES) ?
	                                freopen(currlist[i_open],"r",stdin) :
	                                fopen(currlist[i_open],"r");
	    if (olist[i_open].Io_type.file == NULL)
	    {
	        debug_print("foutput","Left open_close_file()\n");
	        return	NULL;
	    }
	    if (olist[i_open].check_for_blocked_output)
	    {
		IOUTPUT *oput;
	        olist[i_open].output_blocked =
		    check_output(olist[i_open].Io_type.file);
	        oput = index_of_read_file(olist[i_open].Io_type.file,ADD);
		olist[i_open].Io_type = oput->Io_type;
	    }
	    reset_read_file_variables(&olist[i_open]);
	    debug_print("foutput","Left open_close_file()\n");
	    return &olist[i_open].Io_type;
	}
	debug_print("foutput","Left open_close_file()\n");
	return NULL;
}		/*end open_close_file*/


EXPORT	boolean	create_directory(
	const char *dname,
	boolean    all_nodes)
{
	boolean io_node = (is_io_node(pp_mynode())) ? YES : NO;
	boolean status  = FUNCTION_SUCCEEDED;
	char	*path_component, *next_path_component;
	struct stat Buf;
	static	char	*dir = NULL, *partial_path = NULL;
	static	size_t	alloc_len = 0;

	debug_print("directory","Entered create_directory()\n");

	if ( (dname != NULL) && (*dname != '\0') &&
	    ( (all_nodes == YES) || (io_node == YES)) )
	{
	    if (debugging("directory"))
	        (void) printf("Creating directory %s\n",dname);

	    if (alloc_len < strlen(dname))
	    {
	        if (dir != NULL)
	            free(dir);
	        if (partial_path != NULL)
	            free(partial_path);

	        uni_array(&dir,strlen(dname)+1,sizeof(char));
	        uni_array(&partial_path,strlen(dname)+1,sizeof(char));
	        alloc_len = strlen(dname);
	    }

	    (void) strcpy(dir,dname);
	    path_component = strtok(dir,"/");
	    if (dir[0] == '/')
	        (void) sprintf(partial_path,"/%s",path_component);
	    else
	        (void) strcpy(partial_path,path_component);

	    next_path_component = NULL;
	    do
	    {
	        if (next_path_component != NULL)
		{
		    char buff[1024];
		    sprintf(buff,"%s",partial_path);
	            (void) sprintf(partial_path,"%s/%s",buff,next_path_component);
		}
	        if (stat(partial_path,&Buf) != 0)
	        {
	            if (errno == ENOENT)
	            {
	                if (debugging("directory"))
	                {
	                    (void) printf("Creating path component %s\n",
	                                  partial_path);
	                }
	                if (mkdir(partial_path,00755) != 0)
	                {
	                    (void) printf("WARNING in create_directory(), "
	                                  "mkdir() failed for path %s.\n",
	                                  partial_path);
	                    status = FUNCTION_FAILED;
			    break;
	                }
	            }
	            else
	            {
	                (void) printf("WARNING in create_directory(), "
	                              "stat() failed for path %s.\n",
	                              partial_path);
	                status = FUNCTION_FAILED;
	            }
	        }
	        else if (debugging("directory"))
	        {
	            (void) printf("Path component %s already exists\n",
	                          partial_path);
	        }
	    } while ((next_path_component = strtok(NULL,"/")) != NULL);
	}
	/*status = pp_max_status(status); */
	debug_print("directory","Left create_directory()\n");
	return status;
}		/*end create_directory*/

/*
*		set_read_float_size():
*		set_reverse_endian():
*		set_read_endian():
*
*	Provides backward compatibility with older files that did not
*	print the floating point size or endian of the printout.
*/

EXPORT	void	set_read_float_size(
	size_t value)
{
	default_read_float_size = value;
}		/*end set_read_float_size*/

EXPORT	void	set_reverse_endian(
	boolean value)
{
	default_reverse_endian = value;
}		/*end set_reverse_endian*/

EXPORT	void	set_read_endian(
	FT_ENDIAN value)
{
	default_read_endian = value;
}		/*end set_read_endian*/

/*
*			determine_io_type():
*
*	Determines the floating point and endian parameters of an read file.
*/

EXPORT	void	determine_io_type(
	FILE    *file,
	IO_TYPE *io_type)
{
	static OUTPUT *oput = NULL;
	IOUTPUT *ioput;

	if (determine_io_type_blocked)
	    return;

	debug_print("io_type","Entered determine_io_type()\n");
	if ((ioput = index_of_read_file(file,ADD)) == NULL)
	{
	    io_type->file = file;
	    io_type->read_float_size = default_read_float_size;
	    io_type->read_endian     = default_read_endian;
	    io_type->reverse_endian  = default_reverse_endian;
	    io_type->cpu_float_size  = sizeof(double);
	    io_type->cpu_endian = ft_endian_type();
	    if (debugging("foutput"))
	      (void) printf("WARNING in determine_io_type(), "
	                    "index_of_read_file() returns NULL\n");
	    return;
	}
	if (ioput->Io_type.read_endian == FT_UNKNOWN_ENDIAN)
	{
	    determine_io_type_blocked = YES;

	    io_type->file = file;
	    io_type->read_float_size = default_read_float_size;
	    io_type->read_endian     = default_read_endian;
	    io_type->reverse_endian  = default_reverse_endian;
	    io_type->cpu_float_size  = sizeof(double);
	    io_type->cpu_endian = ft_endian_type();
	    oput = save_read_file_variables(file,oput);
	    rewind_read_file(file,NULL);
	    if (next_output_line_containing_string(file,"MACHINE PARAMETERS"))
	    {
	        if (fgetstring(file,"\tByte Ordering            = "))
	        {
	            char line[25];
		    (void) fgets(line,24,file);
		    if (strstr(line,"Big"))
		       io_type->read_endian = FT_BIG_ENDIAN;
		    else if (strstr(line,"Little"))
		        io_type->read_endian = FT_LITTLE_ENDIAN;
		    else
		        io_type->read_endian = FT_UNKNOWN_ENDIAN;
	        }
	        if (fgetstring(file,"\tFloating Point Word Size = "))
	        {
	            (void) fscanf(file,"%lu",&io_type->read_float_size);
	        }
	    }
	    if (((io_type->cpu_endian==FT_BIG_ENDIAN) &&
	         (io_type->read_endian==FT_LITTLE_ENDIAN))
	        ||
	        ((io_type->cpu_endian==FT_LITTLE_ENDIAN) &&
	         (io_type->read_endian==FT_BIG_ENDIAN)))
	        io_type->reverse_endian = YES;
	    else if (((io_type->cpu_endian==FT_BIG_ENDIAN) &&
	              (io_type->read_endian==FT_BIG_ENDIAN))
	        ||
		     ((io_type->cpu_endian==FT_LITTLE_ENDIAN) &&
		      (io_type->read_endian==FT_LITTLE_ENDIAN)) )
	        io_type->reverse_endian = NO;
	    if (io_type == &ioput->Io_type)
	    {
		IO_TYPE Sav_Io_type = *io_type;
	        reset_read_file_variables(oput);
		*io_type = Sav_Io_type;
	    }
	    else
	        reset_read_file_variables(oput);
	    if (debugging("io_type"))
	        (void) printf("Copying io_type to index of read file\n");
	    ioput->Io_type = *io_type;
	    determine_io_type_blocked = NO;
	}
	else
	{
	    if (debugging("io_type"))
	        (void) printf("Copying io_type from index of read file\n");
	    *io_type = ioput->Io_type;
	}

	if (debugging("io_type"))
	    fprint_io_type(stdout,"IO_TYPE from determine_io_type\n",io_type);
	debug_print("io_type","Left determine_io_type()\n");
}		/*end determine_io_type*/

EXPORT	void	fprint_io_type(
	FILE          *file,
	const char    *mesg,
	const IO_TYPE *io_type)
{
        if (mesg != NULL)
	    (void) fprintf(file,"%s",mesg);
	(void) fprintf(file,"IO_TYPE structure %p\n",io_type);
	(void) fprintf(file,"\tfile = %p\n",io_type->file);
	(void) fprintf(file,"\tread_float_size = %llu\n",
	              size_t2uint64_t(io_type->read_float_size));
	(void) fprintf(file,"\tread_endian = %s\n",
	              ft_endian_name(io_type->read_endian));
	(void) fprintf(file,"\treverse_endian = %s\n",
	              y_or_n(io_type->reverse_endian));
	(void) fprintf(file,"\tcpu_float_size = %llu\n",
	              size_t2uint64_t(io_type->cpu_float_size));
	(void) fprintf(file,"\tcpu_endian = %s\n",
	              ft_endian_name(io_type->cpu_endian));
	(void) fprintf(file,"End IO_TYPE structure %p\n",io_type);
}		/*end fprint_io_type*/

EXPORT	uint64_t size_t2uint64_t(
	size_t i)
{
	uint64_t l;
	l = i;
	return l;
}		/*end size_t2uint64_t*/

