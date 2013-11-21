#include "pfdebug.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#include <boost/iterator/iterator_concepts.hpp>
#endif




using namespace std;

pfdebug::pfdebug(int argc, char *argv[])
{
    //this->argc = argc;
    //this->argv = argv;
    //this->desc

    //this->fnm defines in object, which is not so well styled(no idea how to improve)

  
    //set number of threads = 1;
    //gmx_omp_set_num_threads(1);
    omp_set_num_threads(1);
    //mdrunner initialize:
    //in line 953
    this->inputrec      =   new t_inputrec;
    this->mtop          =   new gmx_mtop_t;
    this->state         =   new t_state;

    //inprinciple not a good style, but for convinience of imply libaries
    //initial for cr;
    //runner.c mdrunner() 1071 should not use, but no idea how to repeat.
    this->cr            =   new t_commrec;
    this->cr->mpi_comm_mysim = NULL;
    this->cr->mpi_comm_mygroup = NULL;
    this->cr->nnodes = 1;
    this->cr->sim_nodeid = 0;
    this->cr->nodeid    = this->cr->sim_nodeid;
    
    //copy from cr for -np 1
    this->cr->duty 	=3;
    this->cr->nrank_intranode =1;
    this->cr->nrank_pp_intranode = 1;
    
    //Initial donb_flags, for nonbonded calculation
    this->donb_flags = 0;
    this->donb_flags |= GMX_NONBONDED_DO_SR;		// use the short_nblist
    this->donb_flags |= GMX_NONBONDED_DO_FORCE;
    
    //Initail force container
    snew(this->force_bonded,this->state->natoms);
    snew(this->force_internal, this->state->natoms);
    snew(this->force_external, this->state->natoms);
    snew(this->f_global, this->state->natoms);
    
    //Initial graph
    
    //Initial output flags:
    this->mdof_flogs = 0;
    this->mdof_flogs |= MDOF_X;		//debug use, should cancelled
    this->mdof_flogs |= MDOF_F;		//debug use, result to check
    
    
}

pfdebug::~pfdebug()
{
}

int pfdebug::nfile()
{
    return asize(fnm);
}

void pfdebug::parse_comm(int argc, char* argv[])
{
    // in the domain can we visit the private members.
    parse_common_args(&argc, argv, this->getPCA_Flags(), this->nfile(), this->fnm, asize(this->pa), this->pa,
                      asize(this->desc), this->desc, 0, NULL, &this->oenv);
}




unsigned long pfdebug::getPCA_Flags() const
{
    return PCA_Flags;
}

void pfdebug::setPCA_Flags(unsigned long value)
{
    PCA_Flags = value;
}

void pfdebug::read_top()
{
    read_tpx_state(ftp2fn(efTPX, this->nfile(), this->fnm), this->inputrec, this->state, NULL, this->mtop);
    //related to rerunMD
    this->inputrec->nstlist = 1;
    this->inputrec->nstcalcenergy = 1;
    //nstglobalcomm not really necessary;
    
    //Initial output state, for trajectory only.
    //notice: if active the file, remember to touch the -o file first
    //this->outf = init_mdoutf(this->nfile(), this->fnm, MD_APPENDFILES, this->cr, this->inputrec, this->oenv);
    
    
    // initial atom_to_mol as well;
    //atom_to_mol = new vector<int>;
    //remember, for DNA should be redifined.
    for (int i = 0; i < this->mtop->mols.nr; i++)
    {
      for (int j = this->mtop->mols.index[i]; j < this->mtop->mols.index[i+1]; j++)
      {
	//atom_back current = new atom_back;
	//current.mol_id = i;
	if (i < this->mtop->mols.nr - 1)	//for DNA
	  this->atom_to_mol.push_back(i);
	else 
	  this->atom_to_mol.push_back(i-1);	//for DNA
      }
    }
    
    //initial the residue list:
    for (int i = 0; i < this->mtop->nmoltype; i++)
    {
      t_atoms *atams_temp = &this->mtop->moltype[i].atoms;
      
      // initial resdiues list, for each atom contain the residue number in the molecule
      for (int j = 0; j< atams_temp->nr; j++)
      {
	this->atom_to_res.push_back(atams_temp->atom[j].resind);
      }
      
      //push residue name into string list
      for (int j = 0; j < atams_temp->nres; j++)
      {
	this->residue_list.push_back(atams_temp->resinfo[j].name[0]);
      }
    }
    
}

void pfdebug::init_outf()
{
  //initial as trunc, which means over write every thing
  out_force_internal.open ("force_internal.csv", ofstream::out | ofstream::trunc);
  out_force_external.open ("force_external.csv", ofstream::out | ofstream::trunc);
  out_force_bonded.open ("force_bonded.csv", ofstream::out | ofstream::trunc);
  
  out_force_bonded << "Time(ps)";
  out_force_external << "Time(ps)";
  out_force_internal << "Time(ps)";
  
  int loop_time = this->residue_list.size();
  for (int i = 0; i< loop_time; i++)
  {
    // residue number start from 1, so output start from one.
    out_force_external << ","<< i+1 << "-" << this->residue_list[i] ;
    out_force_internal << ","<< i+1 << "-" << this->residue_list[i] ;
    out_force_bonded << ","<< i+1 << "-" << this->residue_list[i] << ".x," << i+1 << "-"<< this->residue_list[i] << ".y," << i+1 << "-" << this->residue_list[i] << ".z";
  }
  out_force_external << endl;
  out_force_internal << endl;
  out_force_bonded << endl;
  
  out_force_bonded.close();
  out_force_external.close();
  out_force_internal.close();
}


void pfdebug::set_entries()
{
    //set line = ??
    set_state_entries(this->state, this->inputrec, 1);
}



void pfdebug::init_force()
{
    this->fcd           = new t_fcdata;   //this is a pointer, to check if it will work with out knowing the size.

    //runner
    init_disres(this->fplog, this->mtop, this->inputrec, NULL, 0, this->fcd, this->state, false);   //FLAGS = 0

    //runner 1267 
    //read in log file
    gmx_log_open(ftp2fn(efLOG, this->nfile(), this->fnm), this->cr, 1, 0, &this->fplog);	//bMasterOnly =1, bAppendfiles = 0; 
    
    //runner 1371 failed;
    //this->wcycle = wallcycle_init(this->fplog, -1, NULL, 1, 1);

    //runner 1382
    this->nrnb          = new t_nrnb;

    //runner 1400
    this->fr            = new t_forcerec;
    this->fr = mk_forcerec();
    
    //this->fr->hwinfo = hwinfo;
    //set fr->excl_load;
    init_forcerec(this->fplog,this->oenv,this->fr,this->fcd,this->inputrec,this->mtop,this->cr,this->box,
      NULL,NULL,NULL,NULL,"auto",FALSE,-1
    );
    //try to set nthreads == 1;
    this->fr->nthreads = 1;

    //runner 1427
    this->mdatoms = init_mdatoms(this->fplog, this->mtop, this->inputrec->efep != efepNO);

    //runner 1435 vsite
    /*
    this->vsite = init_vsite(this->mtop, NULL, FALSE);

    //runner 1450
    if(this->vsite)
    {
        construct_vsites_mtop(this->vsite, this->mtop,this->state->x);
    }
    */

    //1402
    //fr->hwinf = ??
    //init_forcerec()??
}

void pfdebug::init_fake_md()
{
    //md.c 310
    /*
    init_md(this->fplog,this->cr,this->inputrec,this->oenv,&this->t,&this->t0,this->state->lambda,
            &(this->state->fep_state),this->lam0,
            this->nrnb,this->mtop,&this->upd,
            this->nfile(),this->fnm,&this->outf, &this->mdebin,
            this->force_vir,this->shake_vir,this->mu_tot,&this->bSimAnn,&this->vcm,this->Flags);
            */
    //301
    this->groups = &this->mtop->groups;		//mtop is top global
    
    //md.c 303??
    //init_md()
    
    
    clear_mat(this->total_vir);
    clear_mat(this->pres);
    this->enerd = new gmx_enerdata_t;
    init_enerdata(this->mtop->groups.grps[egcENER].nr, this->inputrec->fepvals->n_lambda,this->enerd);

    //md.c 322
    snew(this->f, this->mtop->natoms);                      //atomic level force
    
    
    //New: made it:
    snew(this->force_bonded, this->mtop->natoms);
    snew(this->force_internal,this->mtop->natoms);
    snew(this->force_external, this->mtop->natoms);
    
    
    
    
    //f = new rvec[mtop->natoms];
    
    //md.c 330
    //copy_df_history()

    //something for energy??
    
    //md.c 342
    //global state
    this->gstat = global_stat_init(this->inputrec);
    
     //md.c 399
    this->top = gmx_mtop_generate_local_top(this->mtop, inputrec);
    
    //md.c 405
    forcerec_set_excl_load(this->fr, this->top, this->cr); //seems global property

    //md.c 408
    this->f_global = this->f;

    //md.c 410  set all atoms to mdatoms, 
    //may change to select atoms => NULL to index
    atoms2md(this->mtop, this->inputrec, 0,NULL, 0, this->mtop->natoms,this->mdatoms);
    if(this->vsite)
    {
      set_vsite_top(this->vsite, this->top, this->mdatoms, this->cr);
    }
}

void pfdebug::md_loop()
{
    //start from md.c 
    //line 732
    this->rerun_fr.natoms = 0;
    
    //md.c 779
    calc_shifts(this->rerun_fr.box, fr->shift_vec);
    
   
    

    //md.c 749
    //readin frame, without checking the natomsnumber
    this->bNotLastFrame = read_first_frame(this->oenv,&this->status,opt2fn("-f",this->nfile(),this->fnm),
                                           &this->rerun_fr, TRX_NEED_X | TRX_READ_V);

    //?? md.c 787 - 799 initialize, prepare the main loop
    this->bFirstStep    = TRUE;
    this->bLastStep     = FALSE;

    this->step  = this->inputrec->init_step;
    this->step_rel = 0;

    //the main loop of MD
    //line 816 to 2178
    //make sure that bNotLastFrame is going to be renewed.
    while(this->bNotLastFrame)
    {
        //cout << this->step << "this is the marker" << endl;
        //819
        wallcycle_start(this->wcycle, ewcSTEP);

        //821 - 842
        this->step      = this->rerun_fr.step;
        this->step_rel  = this->step - this->inputrec->init_step;
        this->t  = this->rerun_fr.time;

        //859-916
            //890 -891 two lines
        copy_mat(this->rerun_fr.box, this->state->box);
	
	//TODO: construct_vsites()

        //931, step for nblist, may not necessary.
        this->bNS = (this->bFirstStep || this->inputrec->nstlist != 0);
        this->bNSList = this->bNS;

        //1028 - 1038 ??

        //1039  cals kinetic

        clear_mat(this->force_vir);

        //1126 init force flag
        //Todo it.
        this->force_flags = 0;

        this->do_me_force();//!!!
	
	//force has been done, then output 
	//md.c 1411
	//this->corretness_test();
	
	
	//write_file
	this->write_file();

	//After write, file, update runtime information
	//this.run_time = gmx_gettime() - (double)runtime->real;

        //this->bNotLastFrame = FALSE;
	//md.c 2088
	this->bNotLastFrame = read_next_frame(this->oenv, this->status, &this->rerun_fr);
    }



}

void pfdebug::corretness_test()
{
  rvec sum_f_bond;
  rvec sum_f_external;
  rvec sum_f_internal;
  clear_rvec(sum_f_bond);
  clear_rvec(sum_f_external);
  clear_rvec(sum_f_internal);
  for (int i = 0; i < this->mtop->natoms; i++)
  {
    rvec_add(this->force_internal[i],sum_f_internal,sum_f_internal);
    rvec_add(this->force_bonded[i],sum_f_bond,sum_f_bond);
    rvec_add(this->force_external[i],sum_f_external,sum_f_external);
  }
  
  for (int i = 0; i < DIM; i++)
  {
    cout << "Dim " << i << "of bond/internal/external summed force is:"<< endl;
    cout << sum_f_bond[i] << endl;
    cout << sum_f_internal[i]<< endl;
    cout << sum_f_external[i]<< endl;
    
  }
}


void pfdebug::write_file()
{
  out_force_internal.open ("force_internal.csv", ofstream::out | ofstream::app);
  out_force_external.open ("force_external.csv", ofstream::out | ofstream::app);
  out_force_bonded.open ("force_bonded.csv", ofstream::out | ofstream::app);
  out_force_bonded << this->t << ",";
  out_force_external << this->t << ",";
  out_force_internal << this->t << ",";
  
  int loop_time = this->residue_list.size();
  int atom_past = 0;
  
  rvec force_bond_res;
  rvec force_internal_res;
  rvec force_external_res;

  
  // output residue based force, for bonded force output the vector force
  for (int i = 0; i < this->mtop->natoms; i++)
  {
    if (i == 0)//start pointer
    {
	atom_past= i;
        clear_rvec(force_bond_res);
	clear_rvec(force_internal_res);
	clear_rvec(force_external_res);
	rvec_add(this->force_bonded[i], force_bond_res, force_bond_res);
	rvec_add(this->force_external[i], force_external_res, force_external_res);
	rvec_add(this->force_internal[i], force_internal_res, force_internal_res);
    }
    else if (atom_to_res[i] == atom_to_res[atom_past] )	//the same residue number
    {
	atom_past = i;
	rvec_add(this->force_bonded[i], force_bond_res, force_bond_res);
	rvec_add(this->force_external[i], force_external_res, force_external_res);
	rvec_add(this->force_internal[i], force_internal_res, force_internal_res);
    }
    else		//to the next residue number
    {
      double f_external = sqrt(force_external_res[0]*force_external_res[0] + force_external_res[1]*force_external_res[1]+ force_external_res[2]*force_external_res[2]);
      double f_internal = sqrt(force_internal_res[0]*force_internal_res[0] + force_internal_res[1]*force_internal_res[1]+ force_internal_res[2]*force_internal_res[2]);
      out_force_external << f_external << ",";
      out_force_internal << f_internal << ",";
      out_force_bonded << force_bond_res[0] << "," <<force_bond_res[1] << "," << force_bond_res[2] << ",";
      
      //restart residue information
      atom_past= i;
        clear_rvec(force_bond_res);
	clear_rvec(force_internal_res);
	clear_rvec(force_external_res);
	rvec_add(this->force_bonded[i], force_bond_res, force_bond_res);
	rvec_add(this->force_external[i], force_external_res, force_external_res);
	rvec_add(this->force_internal[i], force_internal_res, force_internal_res);
    }
    

  }
  
      //after that output the final residue information
      double f_external = sqrt(force_external_res[0]*force_external_res[0] + force_external_res[1]*force_external_res[1]+ force_external_res[2]*force_external_res[2]);
      double f_internal = sqrt(force_internal_res[0]*force_internal_res[0] + force_internal_res[1]*force_internal_res[1]+ force_internal_res[2]*force_internal_res[2]);
      out_force_external << f_external << endl;
      out_force_internal << f_internal << endl;
      out_force_bonded << force_bond_res[0] << "," <<force_bond_res[1] << "," << force_bond_res[2] << endl;
      out_force_external.close();
      out_force_internal.close();
      out_force_bonded.close();
  //debug set:
  //output f:
  //clear_rvecs(this->mtop->natoms,this->f_global);
  
  /*
  //sum all the force by using ugly 
  for (int i = 0; i< this->mtop->natoms; i++)
  {
    rvec_add(this->force_internal[i],this->force_external[i],this->f_global[i]);
    rvec_add(this->f_global[i], this->force_bonded[i], this->f_global[i]);
  }
  */
  
  //equals to md.c 1489, write_traj
  /*
  fwrite_trn(this->outf->fp_trn, this->step, this->t, *this->state->lambda, 
	    this->state->box, this->mtop->natoms, this->state->x, NULL, this->f_global);
	    */
  
  
  //剩下的，自己写，residue based file out put
  clear_rvecs(this->mtop->natoms, this->force_bonded);
  clear_rvecs(this->mtop->natoms, this->force_external);
  clear_rvecs(this->mtop->natoms, this->force_internal);
}


//force.c 175
void pfdebug::do_me_force()
{
    //in force.c
    
    //parameters  used:
  
    //this->fplog 	: 		FILE*
    //this->cr		:		t_commrec*
    //this->inputrec	:		t_inputrec*
    //this->step	:		gmx_large_int_t
    //this->nrnb	:		t_nrnb*
    //this->wcycle	:		gmx_wallcycle_t
    //this->top		:		gmx_localtop_t*		should be define 
    //this->groups	:		gmx_groups_t*
    //this->box		:		matrix
    //this->x		:		rvec
    //this->hist	:		history_t
    //this->f		:		rvec
    //this->vir_force	:		tensor
    //this->mdatoms	:		t_mdatoms
    //this->enerd	:		gmx_enerdata_t
    //this->fcd		:		t_fcdata		bonded force
    //this->lambda	:		real*
    //this->graph	:		t_graph*		What the fuck?
    //this->fr		:		t_forcerec*		force rec parameters
    //this->vsite 	:		gmx_vsite_t*		
    //this->mu_tot 	:		rvec
    //this->t		:		double			should be time
    //this->field	:		FILE*			what the fuck?
    //this->ed		:		gmx_edsam_t		what the fuck?
    //this->bBornRadii	:		gmx_bool		for implicit, useless
    //this->flags	:		int			force flag
    
    //do_force(fplog, cr, ir, step, nrnb, wcycle, top, gropus,
    //		state->box, state->x, &state->hist, 
    //		f, force_vir, mdatoms, enerd, fcd,
    //		state->lambda, graph,
    //		fr, vsite, mu_tot, t, outf->fp_field, ed, bBornRadii,
    //		(bNS ? GMX_FORCE_NS : 0) | force_flags)
  
    //then, for cutoff scheme, choose:
    //do_force_cutsGROUP()
    //parameters used:
    //the same as above, 
  
    //line 1674: ns()!!!!
    //parameters used:
    //this->fp		:	FILE*	
    //this->fr		:	t_forcerec*
    //this->matrix	:	box		??from state??
    //groups:		:	gmx_groups_t*	??from where?
    //top 		:	gmx_localtop_t*	
    //md		:	t_mdatoms*
    //nrnb		:	t_nrnb*
    //bFillGrid		:	gmx_bool
    //bDoLongRangeNS	:	gmx_bool
    
  
    /////////do_force_lowlevel()
  
    // reset free energy components.
    // may not necessary.
    // 187
    /*
    for (int i = 0; i< efptNR; i++)
    {
        this->dvdl_nb[i] = 0;
        this->dvdl_dum[i] = 0;
    }
    */

    
    set_pbc(&this->pbc, this->fr->ePBC, this->rerun_fr.box);
    //this->graph = new t_graph;
    this->graph = mk_graph(this->fplog, &(top->idef), 0, this->mtop->natoms, FALSE, FALSE);
    //265 start calculate nb force.
    
    //donb_flags set here. do_sr | force | energy | do_lr

    //wallcycle_sub_start(wcycle, ewcsNONBONDED);
    
    //sim_util.c line 1542
    //else if (bCalcCGCM)
    calc_cgcm(this->fplog, 0, this->top->cgs.nr, &this->top->cgs, this->rerun_fr.x,this->fr->cg_cm);
    inc_nrnb(this->nrnb, eNR_CGCM, this->mtop->natoms);
    
    //sim_util.c 1675
    //ns()
    this->ns_pf();
    
    //after ns, comes the MAGIC TIME;
    //Sequncial Jobs:
    //1 modify the nblist, only left internal(modified file, use a bool?)
    //2 call do_nonbonded(), force will be internally
    //3 modify the nblist, only left external
    //4 call do_nonbonded(), force will be externally
    
    //Lets go:
    //step 0: do the normal nonbonded:
    //sucessfull!!
    
    /*
    do_nonbonded(this->fr, this->rerun_fr.x, this->force_internal, NULL,this->mdatoms, this->excl,
      &this->enerd->grpp, this->nrnb, this->lambda, this->dvdl_nb, -1, -1, this->donb_flags
    );
    */
    
    //step 1: modify the nblist, only left internal
    // initial nb_list_sr and nb_list_lr; once done, forever done.
    
    this->nb_list_initial();
    
    //step 2 set nblist to internal and call do_nonbonded();
    bool isInternal = true;
    this->ns_modifier(isInternal);
    this->do_the_nonbonded(this->force_internal);
    
    //step 3 and 4 
    
    isInternal = false;
    this->ns_modifier(isInternal);
    this->do_the_nonbonded(this->force_external);    
    

    //step N;
    // destroy the nb_list after used;
    this->nb_list_destroy();
    
    //After the nonbonded, next, we do bonded
    // first, try use the normal calc_bonds function
    

    calc_bonds(this->fplog, this->cr->ms,
	       &this->top->idef, this->rerun_fr.x, &this->state->hist, this->force_bonded, 
	       this->fr, &this->pbc, this->graph,
	       this->enerd, this->nrnb, this->state->lambda, this->mdatoms, this->fcd, NULL, 
	       &this->top->atomtypes,NULL,
	       GMX_FORCE_BONDED | GMX_FORCE_FORCES , false, this->step);

    //504 ewald_LRcorrectiion
    //Above all ,force seems ok, then should do the out put.for the three force categories.
    //No Evald correction is needed.
    
    //out of do me force, go back to md_loop()
    
}

void pfdebug::do_the_nonbonded(rvec* f)
{
  
  do_nonbonded(this->fr, this->rerun_fr.x, f, NULL, this->mdatoms, this->excl,
                 &this->enerd->grpp, this->nrnb, this->lambda, this->dvdl_nb, -1, -1, this->donb_flags);
}

void pfdebug::ns_pf()
{
  // sim_util 1507
  update_forcerec(this->fr, this->state->box);
  
  
  wallcycle_start(this->wcycle, ewcNS);
 
  mk_mshift(this->fplog, this->graph, this->fr->ePBC, this->rerun_fr.box, this->rerun_fr.x);
  //use normal ns function to do the nb search.
  
  
  
  // remember to renew state!!
  //sim_util.c 1674
  ns(this->fplog, this->fr, this->state->box, this->groups,this->top, this->mdatoms,this->cr,this->nrnb,1,0);

  wallcycle_stop(wcycle, ewcNS);
}

void pfdebug::nb_list_initial()
{
  //force.c 83
  this->fr->ns.nblist_initialized = 0;
  int nblength = 10;		// actually waste, but who care
  //this->nb_list_lr.~vector();
  //this->nb_list_sr.~vector();
  t_nblist *nblist_debug;

  // as the container can only be length of 10, never minded;
  for (int k = 0; k < nblength; k++)
  {
    nblist_debug = &this->fr->nblists->nlist_sr[k];
    vector<atom_id> short_list(nblist_debug->nrj, -4);
    vector<atom_id> long_list(nblist_debug->nrj, -4);
      for (int i = 0; i < nblist_debug->nri ; i++) // i is not the atom id
      {
	/*
	cout << i << " Atom in nblist iinr is: " << nblist_debug->iinr[i] ; // this is atom id
	cout << " gid is :" << nblist_debug->gid[i];
	cout << " shift is: " << nblist_debug->shift[i];
	cout << " j range starts at " << nblist_debug->jindex[i];
	cout << " j range ends at " << nblist_debug->jindex[i+1];
	cout << endl;
	*/
	
    
	atom_id i_check = nblist_debug->iinr[i];
	int j_short = nblist_debug->jindex[i];
	int j_long = nblist_debug->jindex[i];
    
	for (int j = nblist_debug->jindex[i]; j< nblist_debug->jindex[i+1] && nblist_debug->jjnr[j] >= 0; j++)
	{
	  int j_check = nblist_debug->jjnr[j];
	  if (atom_to_mol[j_check] == atom_to_mol[i_check])
	  {
	    short_list[j_short++] = j_check;
	    //j_short++;
	  }
	  else
	  {
	    long_list[j_long++] = j_check;
	    //j_long++;
	  }
	}
    
    //debug style
    /*
    cout << "for i atom "<< i << ":" << endl;
    cout << "range is :" << nblist_debug->jindex[i] << " to " << nblist_debug->jindex[i+1] << endl;  
      for (int j = nblist_debug->jindex[i]; j < nblist_debug->jindex[i+1]; j++)
      {
	//cout << endl;
	cout << nblist_debug->jjnr[j] << " "; 
      }
      
      cout << endl;
      
      cout << "short list";
      for (int j = nblist_debug->jindex[i]; j < nblist_debug->jindex[i+1]; j++)
      {
	//cout << endl;
	//cout << jjnr_short[j] << " ";
	cout << short_list[j] << " ";
      }
      cout << endl;
      
      cout << "long list";
      for (int j = nblist_debug->jindex[i]; j < nblist_debug->jindex[i+1]; j++)
      {
	//cout << endl;
	//cout << jjnr_short[j] << " ";
  
	cout << long_list[j] << " ";
      }    
      cout << endl;
      */
      
    }
    
  
    nb_list_lr.push_back(long_list);
    nb_list_sr.push_back(short_list);
  }
}

void pfdebug::nb_list_destroy()
{
  //this->nb_list_lr.~vector();
  //this->nb_list_sr.~vector();
  for (int i = 0; i < 10; i++)
  {
    this->nb_list_lr.pop_back();
    this->nb_list_sr.pop_back();
    this->fr->nblists->nlist_sr[i].jjnr = NULL;// this is important
    //this->nb_list_lr.erase(this->nb_list_lr.begin(), this->nb_list_lr.end());
    //this->nb_list_sr.erase(this->nb_list_sr.begin(), this->nb_list_sr.end());
  }
  //cout << "the size of vector is: " << this->nb_list_lr.size() << endl;
}

    
void pfdebug::ns_modifier(bool isInternal)
{
  /*
   * t_nblist: data structure:
    int igeometry; 	// The type of list(atom, water, etc.)
    int elec; 		// Coulomb loop type index for kernals
    int ielecmod; 	//Coulomb modifier (e.g. switch/shift) what the hell
    int ivdw; 		//VdW loop type index for kernels
    int ivdwmod; 	//VdW modifier (e.g. switch/shift)
    int type; 		//Type of interaction, listed in gmx_nblist_interaction_type
    int nri, maxnri; 	//Current/max number of i particles
    int nrj, maxnrj; 	//Current/max number of j particles
    int maxlen; 	//maxnr of j atoms for a single i atom
    int* iinr; 		//The i-elements
    int* iinr_end; 	// the end atom only with enlist CG;
    int* gid; 		// Index in energy arrays
    int* shift; 	// Shift vector index;
    int* jindex; 	//index in jjnr
    int* jjnr; 		//The j-atom list
    int* jjnr_end; 	//The end atom, only with enltypeCG

    t_excl* excl; 	//Exclustions, only with enltypeCG

    void* kernelptr_vf;
    void* kernelptr_v;
    void* kernelptr_f;
    int simd_padding_width;
   */
  //try to set nlist_lr = nlist_sr
  //failed: invalid array assignment.
  //this->fr->nblists->nlist_lr = this->fr->nblists->nlist_sr;
  
  //modified result after done the normal ns();
  //ns_debug();
  // Because nblist_sr only have length of 10, so loop for 10 is ok 
  
  int nblength = 10;
  if(isInternal)
  {
    // do the short
    for (int i = 0; i < nblength; i++)
    {
      this->fr->nblists->nlist_sr[i].jjnr = &this->nb_list_sr[i].front();
    }
  }
  else
  {
    // do the long
    for (int i = 0; i < nblength; i++)
    {
      this->fr->nblists->nlist_sr[i].jjnr = &this->nb_list_lr[i].front();
    }
    
  }
}

void pfdebug::moltype_debug()
{
  //Datatype: 			gmx_mtop_t:	The Global topology format of *.tpr
  
  //Including:
  //char			**names
  //gmx_ffparams_t		ffparams:	force field
  //int 			nmoltype	should be the number of molecules
  //gmx_moltype_t:		*moltype	
  //int 			nmolblock:	
  //gmx_molblock_t:		*molblock;
  //int 			natoms;
  //int				maxres_renum;
  //int				maxresnr;	The maximum residue number in moltype
  //t_atomtypes			atomtypes;	Atomtype properties
  //t_block			mols		the molecules
  //gmx_gruops_t		groups;		Charge groups?
  //t_symtab			symtab;		the symbol table
  
  
  //names:
  cout << *this->mtop->name << endl; //easy, cout, Good Job!
  
  //ffparams:		who care
  
  //nmoltype:		
  //how many kind of molecules in *tpr, (test case:2 plus water and salt is 4)
  
  //moltype:	should have the information for each nmoltype:
  /**/
  for (int i = 0; i< this->mtop->nmoltype; i++)
  {
    gmx_moltype_t *temp_mol_i = &this->mtop->moltype[i];
    
    cout << *temp_mol_i->name << " has ";
    cout << temp_mol_i->atoms.nr << " Atoms: " << endl;}
    /*for (int j = 0; j < temp_mol_i->atoms.nr; j++)
    {
      cout << temp_mol_i->atoms.atom[j].atomnumber << " in Res "; 		// atomnumber is the number in the element table!!
      cout << temp_mol_i->atoms.atom[j].resind << endl;
    }
    cout << endl;
    
    cout << temp_mol_i->atoms.nres << " resdiues " << endl;
    for (int j = 0; j< temp_mol_i->atoms.nres; i++)
    {
      //cout << this->mtop->moltype.atoms->resinfo[j].name << " with: nr = " << this->mtop->moltype.atoms->resinfo[j].nr << endl;
    }
    cout << endl;
  }
  */
  
  //nmoltype and nmolblock, no idea
  
  //natoms: number of atoms in the tpr
  
  //maxres_renum and maxresnr no idea, seems the max res number in a mol
  
  //atomtypes:	atom information, normally nr = 19
  
  //mols:
  //index list for each molecules:
  /*
  for (int i = 0; i< this->mtop->mols.nr; i++)
  {
    cout << this->mtop->mols.index[i] << endl;
  }
 */ 
  
  /*
  // mtop->mols contain mols ids according to molecules
  
  cout << this->mtop->mols.nr << endl;
  for (int i = 0; i< 2000; i++)
  {
    cout << this->mtop->mols.index[i] << endl;
  }
  
  
  
  int nt = 1; //moltype
  for(int i = 0; i< this->mtop->moltype[nt].atoms.nr; i++)
  {
    cout << this->mtop->moltype[nt].atoms.atom[i].resind << " atomnumber: ";
    cout << this->mtop->moltype[nt].atoms.atom[i].atomnumber << " particle type";
    cout << this->mtop->moltype[nt].atoms.atom[i].type;
    cout << endl;
    //cout << this->mtop->moltype[nt].nmol;
    //cout << this->mtop->moltype[nt]
    
    //so atomnumber is local in res, not globally
    //cout << "atom id:" << this->mtop->moltype[nt].atoms.atom[i].resind << endl;
  }
  */
  cout << endl;
}


//not use in real simulation
void pfdebug::ns_debug()
{
  //check out what is done and what is not;
  
  t_nblist	*nblist_debug;
  
  
  
  //replace the nlist_lr, somehow
  this->fr->nblists->nlist_lr[0] = this->fr->nblists->nlist_sr[2];
  
  nblist_debug 	= &this->fr->nblists->nlist_sr[2];
  
  int nlist_i = 0;
  //test redefine jjnr;
  
  vector<atom_id> short_list(nblist_debug->nrj, -4);
  atom_id * jjnr_short = &short_list.front(); 
  //atom_id *jjnr_short = new atom_id[nblist_debug->nrj](-4);
  vector<atom_id> long_list(nblist_debug->nrj, -4);
  atom_id * jjnr_long = &long_list.front(); 
  
  this->fr->nblists->nlist_lr[nlist_i].jjnr = jjnr_long;
  
  for (int i = 0; i < nblist_debug->nri ; i++) // i is not the atom id
  {
    /*
    cout << i << " Atom in nblist iinr is: " << nblist_debug->iinr[i] ; // this is atom id
    cout << " gid is :" << nblist_debug->gid[i];
    cout << " shift is: " << nblist_debug->shift[i];
    cout << " j range starts at " << nblist_debug->jindex[i];
    cout << " j range ends at " << nblist_debug->jindex[i+1];
    cout << endl;
    */
    
    atom_id i_check = nblist_debug->iinr[i];
    int j_short = nblist_debug->jindex[i];
    int j_long = nblist_debug->jindex[i];
    
    for (int j = nblist_debug->jindex[i]; j< nblist_debug->jindex[i+1] && nblist_debug->jjnr[j] >= 0; j++)
    {
      int j_check = nblist_debug->jjnr[j];
      if (atom_to_mol[j_check] == atom_to_mol[i_check])
      {
	jjnr_short[j_short++] = j_check;
	//j_short++;
      }
      else
      {
	jjnr_long[j_long++] = j_check;
	//j_long++;
      }
    }
    
    /*
    for (int j = nblist_debug->jindex[i]; j < nblist_debug->jindex[i+1]; j++)
    {
      //cout << endl;
      cout << nblist_debug->jjnr[j] << " "; 
    }
    
    cout << endl;
    for (int j = nblist_debug->jindex[i]; j < nblist_debug->jindex[i+1]; j++)
    {
      //cout << endl;
      //cout << jjnr_short[j] << " ";
      cout << this->fr->nblists->nlist_lr[nlist_i].jjnr[j] << " ";
    }
    */
    
    cout << endl;
  }
  
  cout << nlist_i;
  //cout << this->fr->nblists->nlist_sr
}








//test case
void pfdebug::readin_sec()
{
    std::cout << this->mtop->mols.nr << std::endl;          //output numbers of mols
    cout << this->mtop->natoms << endl;                     //output numbers of atoms
    cout << this->state->natoms << endl;
    cout << this->rerun_fr.natoms << endl;
    //cout << can_use_allvsall(this->inputrec, true, NULL, fplog) << endl;
    //cout << sizeof f << endl;
}







//main function
int main(int argc, char *argv[])
{
    //Sec.1: initialize all the parameters.
    pfdebug *test = new pfdebug(argc, argv);

    //test->setPCA_Flags(PCA_CAN_SET_DEFFNM | 0);
    test->parse_comm(argc,argv);


    //Sec.2: then go inside mdrunner(), complicate, to be abstracted.
    //the debug of trr outpu is set inside.
    test->read_top();       //readin topology to mtop, state, inputrec
    
    
    //Sec.2.1 init output file here.
     test->init_outf();
    
    //test->moltype_debug();
    //test->set_entries();    //??set entries??
    test->init_force();

    //Sec.3 do_md()
    test->init_fake_md();

    //Sec.3.1 md.c 734 Loop over MD steps
    //Inprinciple should be 
    test->md_loop();



    //Sec N: How to dump out file?
    //try in moltype_debug
    //test->moltype_debug();




    //
    cout << "hello gmx!" << endl;
    cout<< test->nfile() << endl;
    test->readin_sec();
    return 0;
}












