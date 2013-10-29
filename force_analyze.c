#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

#include "gromacs/legacyheaders/checkpoint.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/filenm.h"
#include "gromacs/legacyheaders/gmx_fatal.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/readinp.h"
#include "gromacs/legacyheaders/statutil.h"
#include "gromacs/legacyheaders/typedefs.h"

#include "gromacs/commandline/cmdlinemodulemanager.h"

#include "force_analyze.h"

int force_analyze(int argc, char *argv[])
{
  const char *desc[] = {
  "force analyze tool:",
  "is to analyze force."
  };
  t_commrec *cr;    //cr
  t_filenm  fnm[] = {
    { efTPX, NULL, NULL, ffREAD },
  }
}

