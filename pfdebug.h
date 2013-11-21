/*
    Author: Chen Chen
    E-mail: maverickcc@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef PFDEBUG_H
#define PFDEBUG_H

#include <stdio.h>
#include <vector>
#include <string.h>
#include <iostream>

#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/filenm.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/statutil.h"     //parse_common_args
#include "gromacs/legacyheaders/typedefs.h"

#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "filenm.h"

//for mdrunner
#include "tpxio.h"
#include "nrnb.h"
#include "force.h"
#include "disre.h"
#include "gmx_wallcycle.h"
#include "mdatoms.h"

//for do_md
#include "smalloc.h"
#include "mdatoms.h"
#include "mtop_util.h"

//for ns()
#include "force.h"
#include "chargegroup.h"
#include "ns.h"
#include "nrnb.h"

//for nonbond
#include "nonbonded.h"

//for bonded
#include "bondf.h"

//for output
#include "trnio.h"
#include "sim_util.h"

//only use to kill the omp thread
#include "gmx_omp.h"
#include "omp.h"
#include "gmx_omp_nthreads.h"

//for output file
#include <fstream>



using namespace std;

class pfdebug
{

public:
    pfdebug(int argc, char *argv[]);
    virtual ~pfdebug();


    //get NFILE
    int nfile();        //replace #NFILE

    //parse_commen_args
    void parse_comm(int argc, char *argv[]);

    //PCA_Flags get and set
    unsigned long getPCA_Flags() const;
    void setPCA_Flags(unsigned long value);

    //inside mdrunner
    //readin top file to inputrec, state, and mtop;
    void read_top();
    
    //inite output format
    void init_outf();

    //set state entries:    1122
    void set_entries();

    // open fplog??
    void open_log();

    //initial fcd:          1188
    void init_force();


    //for do_md
    void init_fake_md();

    //loop over MD steps
    void md_loop();

    //do_force
    void do_me_force();
    
    //write_file
    void write_file();

    //ns() inside do_me_force()
    //the same api as original ns();
    void ns_pf();
    
    //initial nb_list_lr and nb_list_sr
    void nb_list_initial();
    
    //destroy nb_list_lr and nb_list_sr
    void nb_list_destroy();

    //seperate internal and exteranl
    void ns_modifier(bool isInternal);
    
    //do no_bonded
    //for f_internal or for f_external;
    void do_the_nonbonded(rvec *f);

    //debug ns
    void ns_debug();

    //try to understand the mtop->moltype
    //nmoltype: number of molecules
    void moltype_debug();

    //test case
    void readin_sec();
    
    //corretness test
    //sum of force_internal should be zero.
    //sum of force_external should be zero
    //som of force_bonded shoudl be zero
    void corretness_test();



private:

    //test cases
    //const char *test = "abcde";
    //const int testint[5] = { 1, 2, 3, 4 ,2};

    //int argc;
    //char **argv;

    //for description
    const char      *desc[2] = {
        "what ever",
        "test it"
    };



    //may not necessary, use for paral.
    t_commrec       *cr;

    //fnm input, include file type.
    t_filenm        fnm[9] = {
        { efTPX,        NULL,      NULL,       		ffREAD },
        { efTRN,        "-o",      NULL,       		ffWRITE },
        { efCPT,        "-cpi",    NULL,       		ffOPTRD },
        { efSTO,        "-c",      "confout",  		ffWRITE },
        { efTRX,        "-f",      NULL,       		ffREAD },
        { efLOG,        "-g",      "md",       		ffWRITE},
	{ efXVG,	"-bf",	   "bound_force", 	ffOPTWR},
	{ efXVG,	"-if",	   "internal_force",	ffWRITE},
	{ efXVG,	"-xf",     "external_force", 	ffWRITE}
	
    };

    // Command line options
    gmx_bool      bPartDec      = FALSE;
    gmx_bool      bDDBondCheck  = TRUE;
    gmx_bool	  bDDBondComm   = TRUE;
    gmx_bool      bTunePME      = TRUE;
    gmx_bool      bTestVerlet   = FALSE;
    gmx_bool      bVerbose      = FALSE;
    gmx_bool      bCompact      = TRUE;
    gmx_bool      bSepPot       = FALSE;
    gmx_bool      bRerunVSite   = FALSE;        //might useful
    gmx_bool      bIonize       = FALSE;
    gmx_bool      bConfout      = TRUE;
    gmx_bool      bReproducible = FALSE;

    int           npme          = -1;
    int           nmultisim     = 0;
    int           nstglobalcomm = -1;
    int           repl_ex_nst   = 0;
    int           repl_ex_seed  = -1;
    int           repl_ex_nex   = 0;
    int           nstepout      = 100;
    int           resetstep     = -1;
    gmx_large_int_t nsteps      = -2; /* the value -2 means that the mdp option will be used */

    rvec          realddxyz          = {0, 0, 0};
    const char   *ddno_opt[ddnoNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };

    const char   *dddlb_opt[5] =
    { NULL, "auto", "no", "yes", NULL };

    const char   *thread_aff_opt[threadaffNR+1] =
    { NULL, "auto", "on", "off", NULL };

    const char   *nbpu_opt[6] =
    { NULL, "auto", "cpu", "gpu", "gpu_cpu", NULL };
    real          rdd                   = 0.0, rconstr = 0.0, dlb_scale = 0.8, pforce = -1;
    char         *ddcsx                 = NULL, *ddcsy = NULL, *ddcsz = NULL;
    real          cpt_period            = 15.0, max_hours = -1;
    gmx_bool      bAppendFiles          = TRUE;
    gmx_bool      bKeepAndNumCPT        = FALSE;
    gmx_bool      bResetCountersHalfWay = FALSE;

    //output setup
    output_env_t  oenv                  = NULL;
    const char   *deviceOptions         = "";

    gmx_hw_opt_t  hw_opt = {0, 0, 0, 0, threadaffSEL, 0, 0, NULL};

    // parameters
    t_pargs       pa[1] = {
        {   "-v",       FALSE, etBOOL, {&bVerbose},
            "Be loud and noisy"
        }
    };

    unsigned long Flags;        //Flags for simulation, may not use for analyze.
    unsigned long PCA_Flags;
    ivec          ddxyz;
    int           dd_node_order;
    gmx_bool      bAddPart;
    FILE         *fplog, *fpmulti;
    int           sim_part, sim_part_fn;
    const char   *part_suffix = ".part";
    char          suffix[STRLEN];
    int           rc;
    char        **multidir = NULL;  //may not use for analyze.

    //cr = init_par();

    //PCA_Flags = (PCA_CAN_SET_DEFFNM | (MASTER(cr) ? 0 : PCA_QUIET));
    //others
    FILE            *fp;
    char            *grpname;
    const char      *fdist;
    int             gnx;
    atom_id         *index;
    char            title[STRLEN];
    //t_topology      top;
    int             ePBC = -1;
    rvec            *x;
    matrix          box;

    //parameters in mdrunner:
    t_inputrec      *inputrec;              // all the parameters from tpr file
    gmx_mtop_t      *mtop;                  // global topologies imformation
    gmx_localtop_t  *top;		    // local topology imformation!!!! not yet
    t_state         *state;                 // ??

    gmx_wallcycle_t wcycle;                 //wall molecules.???
    t_mdatoms       *mdatoms;               //?? mass, charge, and other information of all atoms

    gmx_vsite_t     *vsite;                 //vsite


    t_fcdata        *fcd;                   //data type for bonded force.
    t_nrnb          *nrnb;                  //data type for non bonded force(nb list and ele)

    t_forcerec      *fr;                    //should be force rec, the frame information is in rerun_fr

    //parameters for do_md

    gmx_mdoutf_t    *outf;
    gmx_update_t    upd;
    gmx_bool        bSimAnn;
    t_mdebin        *mdebin;
    t_vcm           *vcm;
    double          run_time;               //??
    tensor          force_vir, shake_vir, total_vir, tmp_vir, pres;
    rvec            mu_tot;

    rvec            *f_global;
    rvec            *x_xtc;
    gmx_enerdata_t  *enerd;
    rvec            *f;
    atom_id         *grpindex;
    gmx_groups_t	*groups;
    gmx_global_stat_t	gstat;

    //for main analyse loop
    t_trxframe      rerun_fr;               //for analyze of frame
    t_trxstatus     *status;
    gmx_bool        bNotLastFrame;

    gmx_bool        bFirstStep;
    gmx_bool        bLastStep;

    gmx_bool        bNS;
    gmx_bool        bNSList;

    double          t;
    double          t0;
    double          lam0[efptNR];

    int             force_flags;


    gmx_large_int_t step;
    gmx_large_int_t step_rel;

    //for do_me_force as force.c

    real        dvdl_dum[efptNR], dvdl, dvdl_nb[efptNR];

    //int         donb_flags;
    rvec        *f_longrange;
    t_blocka    *excl;
    real        *lambda;


    //in do_nonbonded

    //NEW:
    
    /*
    typedef struct {
      int mol_id; // atom in mol_id
      //int res_id;
    } atom_back; // index of this type is atom id;
    */
    
    //atoms and res list:
    vector <int> atom_to_mol;
    vector <int> atom_to_res;
    vector <string> residue_list;
    
    //nblist container
    vector <vector <atom_id> >	nb_list_sr;
    vector <vector <atom_id> >  nb_list_lr;
    
    //bonded force
    rvec	*force_bonded;
    t_pbc 	pbc;
    t_graph 	*graph;

    //nonbonded force
    rvec	*force_internal;
    rvec	*force_external;
    int		donb_flags;		// force flag for nb
    
    //for output
    int 	mdof_flogs;
    
    //output the force information
    ofstream	out_force_bonded;
    ofstream	out_force_internal;
    ofstream	out_force_external;
    
};

#endif // PFDEBUG_H
