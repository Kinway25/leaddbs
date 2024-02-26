/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__Fzap
#define _nrn_initial _nrn_initial__Fzap
#define nrn_cur _nrn_cur__Fzap
#define _nrn_current _nrn_current__Fzap
#define nrn_jacob _nrn_jacob__Fzap
#define nrn_state _nrn_state__Fzap
#define _net_receive _net_receive__Fzap 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define del _p[0]
#define del_columnindex 0
#define dur _p[1]
#define dur_columnindex 1
#define f0 _p[2]
#define f0_columnindex 2
#define f1 _p[3]
#define f1_columnindex 3
#define amp _p[4]
#define amp_columnindex 4
#define f _p[5]
#define f_columnindex 5
#define on _p[6]
#define on_columnindex 6
#define _tsav _p[7]
#define _tsav_columnindex 7
#define _nd_area  *_ppvar[0]._pval
#define x	*_ppvar[2]._pval
#define _p_x	_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  2;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "del", "ms",
 "dur", "ms",
 "f0", "1/s",
 "f1", "1/s",
 "amp", "1",
 "f", "1/s",
 "x", "1",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void _ba1(Node*_nd, double* _pp, Datum* _ppd, Datum* _thread, NrnThread* _nt) ;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Fzap",
 "del",
 "dur",
 "f0",
 "f1",
 "amp",
 0,
 "f",
 0,
 0,
 "x",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 8, _prop);
 	/*initialize range parameters*/
 	del = 0;
 	dur = 0;
 	f0 = 0;
 	f1 = 0;
 	amp = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 8;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[3]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _fzap_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 8, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
  hoc_register_dparam_semantics(_mechtype, 3, "netsend");
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_reg_ba(_mechtype, _ba1, 11);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Fzap /home/netstim/leaddbs_git/ext_libs/OSS-DBS/Axon_Processing/Axon_files/McNeal1976/fzap.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define PI _nrnunit_PI[_nrnunit_use_legacy_]
static double _nrnunit_PI[2] = {0x1.921fb54442d18p+1, 3.14159}; /* 3.14159265358979312 */
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 /* BEFORE BREAKPOINT */
 static void _ba1(Node*_nd, double* _pp, Datum* _ppd, Datum* _thread, NrnThread* _nt)  {
    _p = _pp; _ppvar = _ppd;
  v = NODEV(_nd);
 if ( on  == 0.0 ) {
     f = 0.0 ;
     x = 0.0 ;
     }
   else {
     f = f0 + ( f1 - f0 ) * ( t - del ) / dur ;
     x = amp * sin ( 2.0 * PI * ( t - del ) * ( f0 + ( f1 - f0 ) * ( t - del ) / ( 2.0 * dur ) ) * ( 0.001 ) ) ;
     }
   }
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 1.0 ) {
     if ( on  == 0.0 ) {
       on = 1.0 ;
       net_send ( _tqitem, _args, _pnt, t +  dur , 1.0 ) ;
       }
     else {
       on = 0.0 ;
       }
     }
   } }

static void initmodel() {
  int _i; double _save;_ninits++;
{
 {
   f = 0.0 ;
   x = 0.0 ;
   on = 0.0 ;
   if ( del < 0.0 ) {
     del = 0.0 ;
     }
   if ( dur < 0.0 ) {
     dur = 0.0 ;
     }
   if ( f0 <= 0.0 ) {
     f0 = 0.0 ;
     }
   if ( f1 <= 0.0 ) {
     f1 = 0.0 ;
     }
   if ( dur > 0.0 ) {
     net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  del , 1.0 ) ;
     }
   }

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/netstim/leaddbs_git/ext_libs/OSS-DBS/Axon_Processing/Axon_files/McNeal1976/fzap.mod";
static const char* nmodl_file_text = 
  ": $Id: fzap.mod,v 1.5 2015/05/28 02:51:15 ted Exp ted $\n"
  "\n"
  "COMMENT\n"
  "fzap.mod\n"
  "\n"
  "A bogus point process that contains the variable x, \n"
  "which oscillates starting at t = del >= 0.\n"
  "The frequency f of the oscillation increases linearly with time\n"
  "from f0 at t == del to f1 at t == del + dur, \n"
  "where both del and dur are > 0.\n"
  "\n"
  "fzap uses the event delivery system to ensure compatibility with adaptive integration.\n"
  "\n"
  "=================\n"
  "NOTES AND CAVEATS\n"
  "=================\n"
  "\n"
  "1.  If x were a RANGE variable, an assignment statement would \n"
  "have to be inserted into proc advance() in order for the \n"
  "value of x to be used by other mechanisms--e.g.\n"
  "proc advance() {\n"
  "  is_xtra = Fzap[0].x\n"
  "  fadvance()\n"
  "}\n"
  "However, that would be incompatible with adaptive integration.\n"
  "To eliminate the need for such an assignment statement, x is a \n"
  "POINTER.  This preserves compatibility with adaptive integration.\n"
  "\n"
  "2.  On every fadvance, the statements that evaluate Fzap's x \n"
  "should be executed before the statements in any client mechanism \n"
  "that relies on the value of Fzap's x.  To that end, the value of \n"
  "x is computed in a BEFORE BREAKPOINT block, which will take care\n"
  "of any client mechanism that uses Fzap's x in a BREAKPOINT block.\n"
  "\n"
  "However, some client mechanisms may have their own \n"
  "BEFORE BREAKPOINT blocks that need the value of Fzap's x.  \n"
  "xtra is such a mechanism.  In this situation, care is required \n"
  "to ensure that the statements in Fzap's BEFORE BREAKPOINT block\n"
  "are executed first.  This can be done by compiling the mod file \n"
  "that defines Fzap _before_ the client mechanism's mod file.\n"
  "\n"
  "There are two ways to make this happen:\n"
  "A.  Invoke nrnivmodl with a command line that presents the file \n"
  "names in the desired sequence.  UNIX/Linux users may be quite \n"
  "comfortable with this.\n"
  "B.  Choose mod file names so that Fzap's mod file appears before \n"
  "the name of any client mod files in an alphabetical listing.\n"
  "For the example of Fzap and xtra, the file names fzap.mod and \n"
  "xtra.mod would be quite suitable.  This is more convenient for \n"
  "users of all operating systems, but especially MSWin and OS X, \n"
  "whose users are accustomed to compiling all mod files in a \n"
  "directory with mknrndll or \"drag and drop,\" respectively.\n"
  "\n"
  "12/11/2008 NTC\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "  POINT_PROCESS Fzap\n"
  "  RANGE del, dur, f0, f1, amp, f\n"
  "  POINTER x\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "  PI = (pi) (1)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  del (ms)\n"
  "  dur (ms)\n"
  "  f0 (1/s)  : frequency is in Hz\n"
  "  f1 (1/s)\n"
  "  amp (1)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  f (1/s)\n"
  "  x (1)\n"
  "  on (1)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "  f = 0\n"
  "  x = 0\n"
  "  on = 0\n"
  "\n"
  "  if (del<0) { del=0 }\n"
  "  if (dur<0) { dur=0 }\n"
  "  if (f0<=0) { f0=0 (1/s) }\n"
  "  if (f1<=0) { f1=0 (1/s) }\n"
  "\n"
  "  : do nothing if dur == 0\n"
  "  if (dur>0) {\n"
  "    net_send(del, 1)  : to turn it on and start frequency ramp\n"
  "  }\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "The angular velocity in radians/sec is w = 2*PI*f, \n"
  "where f is the instantaneous frequency in Hz.\n"
  "\n"
  "Assume for the moment that the frequency ramp starts at t = 0.\n"
  "f = f0 + (f1 - f0)*t/dur\n"
  "\n"
  "Then the angular displacement is\n"
  "theta = 2*PI * ( f0*t + (f1 - f0)*(t^2)/(2*dur) ) \n"
  "      = 2*PI * t * (f0 + (f1 - f0)*t/(2*dur))\n"
  "But the ramp starts at t = del, so just substitute t-del for every occurrence of t\n"
  "in the formula for theta.\n"
  "ENDCOMMENT\n"
  "\n"
  "BEFORE BREAKPOINT {\n"
  "  if (on==0) {\n"
  "    f = 0\n"
  "    x = 0\n"
  "  } else {\n"
  "    f = f0 + (f1 - f0)*(t-del)/dur\n"
  "    x = amp * sin( 2*PI * (t-del) * (f0 + (f1 - f0)*(t-del)/(2*dur)) * (0.001) )\n"
  "  }\n"
  "}\n"
  "\n"
  "NET_RECEIVE (w) {\n"
  "  : respond only to self-events with flag > 0\n"
  "  if (flag == 1) {\n"
  "    if (on==0) {\n"
  "      on = 1  : turn it on\n"
  "      net_send(dur, 1)  : to stop frequency ramp, freezing frequency at f1\n"
  "    } else {\n"
  "      on = 0  : turn it off\n"
  "    }\n"
  "  }\n"
  "}\n"
  ;
#endif
