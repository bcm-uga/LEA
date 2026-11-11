/* Instrument actual ANCESTRYMAP code: every explicit R boundary must have
 * released owned allocations/streams. Caller-owned buffers remain borrowed. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <setjmp.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
static void *pointers[64];
static FILE *streams[8], *read_error;
static size_t operations, fail_at, reads, fail_read;
static int caught;
static char diagnostic[2048];
static jmp_buf escape;
static int allocation_slot(void *p) { for(int i=0;i<64;i++)if(pointers[i]==p)return i;abort(); }
static void clean(void) {
    for(int i=0;i<64;i++)assert(!pointers[i]);
    for(int i=0;i<8;i++)assert(!streams[i]);
}
static int fault(void) { return ++operations == fail_at; }
static void *tracked_realloc(void *p,size_t n) {
    if(fault())return NULL;
    int i=allocation_slot(p);void *replacement=realloc(p,n);
    if(replacement)pointers[i]=replacement;return replacement;
}
static void *tracked_calloc(size_t n,size_t width) {
    if(fault())return NULL;void *p=calloc(n,width);
    if(p)pointers[allocation_slot(NULL)]=p;return p;
}
static void tracked_free(void *p) {if(p)pointers[allocation_slot(p)]=NULL;free(p);}
static FILE *tracked_fopen(const char *path,const char *mode) {
    if(fault())return NULL;FILE *f=fopen(path,mode);
    if(f){int i=0;while(i<8&&streams[i])i++;assert(i<8);streams[i]=f;}return f;
}
static int tracked_fclose(FILE *f) {
    int failure=fault(),i=0;while(i<8&&streams[i]!=f)i++;assert(i<8);
    streams[i]=NULL;if(read_error==f)read_error=NULL;
    int result=fclose(f);return failure?EOF:result;
}
static int tracked_fseek(FILE *f,long offset,int origin) {return fault()?-1:fseek(f,offset,origin);}
static int tracked_fgetc(FILE *f) {if(++reads==fail_read){read_error=f;return EOF;}return fgetc(f);}
static int tracked_ferror(FILE *f) {return read_error==f||ferror(f);}
static int tracked_fprintf(FILE *f,const char *fmt,...) {
    if(fault())return -1;va_list args;va_start(args,fmt);int result=vfprintf(f,fmt,args);va_end(args);return result;
}
static int tracked_fputc(int c,FILE *f) {return fault()?EOF:fputc(c,f);}
#define realloc tracked_realloc
#define calloc tracked_calloc
#define free tracked_free
#define fopen tracked_fopen
#define fclose tracked_fclose
#define fseek tracked_fseek
#define fgetc tracked_fgetc
#define ferror tracked_ferror
#define fprintf tracked_fprintf
#define fputc tracked_fputc
#include "convert/ancestrymap.c"
#undef realloc
#undef calloc
#undef free
#undef fopen
#undef fclose
#undef fseek
#undef fgetc
#undef ferror
#undef fprintf
#undef fputc
void Rprintf(const char *format,...) { (void)format;clean(); }
_Noreturn void error(const char *format,...) {
    assert(format&&*format);clean();
    va_list args;va_start(args,format);vsnprintf(diagnostic,sizeof(diagnostic),format,args);va_end(args);
    caught=1;longjmp(escape,1);
}
static void reset(size_t target,size_t read_target) {
    operations=reads=0;fail_at=target;fail_read=read_target;caught=0;diagnostic[0]=0;
}
static void run(char *input,char *output,int lfmm,size_t target,size_t read_target,int expected) {
    int n=0,m=0;reset(target,read_target);
    if(!setjmp(escape)) {
        if(lfmm)ancestrymap2lfmm(input,output,&n,&m);else ancestrymap2geno(input,output,&n,&m);
        assert(n==2&&m==2);
    }
    clean();assert(caught==expected);
    if(!expected) {
        char text[128];FILE *f=fopen(output,"r");assert(f);
        size_t nread=fread(text,1,sizeof(text)-1,f);text[nread]=0;assert(!fclose(f));
        assert(!strcmp(text,lfmm?"0 9\n3 -1\n":"03\n9-1\n"));
    }
}
int main(int argc,char **argv) {
    assert(argc==3);
    size_t injected=0;
    for(int lfmm=0;lfmm<2;lfmm++) {
        run(argv[1],argv[2],lfmm,0,0,0);
        size_t total=operations, total_reads=reads;
        for(size_t i=1;i<=total;i++){run(argv[1],argv[2],lfmm,i,0,1);injected++;}
        size_t positions[]={1,2,50,total_reads/2-1,total_reads/2+5,total_reads-1,total_reads};
        for(size_t i=0;i<sizeof(positions)/sizeof(positions[0]);i++) {
            run(argv[1],argv[2],lfmm,0,positions[i],1);injected++;
        }
        run(argv[1],argv[2],lfmm,0,0,0);
    }
    reset(0,0);assert(nb_ind_ancestrymap(argv[1])==2);clean();
    int values[4]={0};read_ancestrymap(argv[1],2,2,values);clean();
    assert(values[0]==0&&values[1]==9&&values[2]==3&&values[3]==-1);
    /* A successful final read must not be discarded when dimensions are full;
     * borrowed stack storage must never be freed at the error boundary. */
    int borrowed=123;reset(0,0);
    if(!setjmp(escape))read_ancestrymap(argv[1],1,1,&borrowed);
    assert(caught&&strstr(diagnostic,"unexpected extra records"));clean();
    am_context ctx={0};ctx.input_name="integer-boundary";
    size_t count=(size_t)INT_MAX-1;
    assert(am_increment(&ctx,&count)&&count==(size_t)INT_MAX);
    assert(!am_increment(&ctx,&count)&&count==(size_t)INT_MAX);
    ctx.N=INT_MAX;ctx.M=2;assert(!am_matrix_size(&ctx));
    ctx.N=0;ctx.M=1;assert(!am_matrix_size(&ctx));
    ctx.N=1;ctx.M=1;assert(am_matrix_size(&ctx)&&ctx.cells==1);
    printf("%zu injected allocation/open/seek/read/write/close failures passed; "
           "zero live resources at R boundaries; borrowed-buffer, extra-record and integer-limit checks passed.\n",injected);
    return 0;
}
