#!/bin/bash

ncpu=6
# ncpu=$SLURM_NTASKS

QEBIN='../../bin'
BIN_DIR=$QEBIN
EXEC="mpirun"
ECHO=echo

lrun_scf=true
lrun_rmta=true

lsave_tmp_dir=false
lsave_tmp_dir_tar=false

PREFIX='nb-bcc'
IBRAV="-3"
ALAT=6.2383834
NAT=1
NTYP=1
ECUT=100
CHARGE=0
PRESS=0
OCCUPATIONS="smearing"
SMEARING='mp'
DEGAUSS=0.01
PSEUDO_DIR='pseudo'
TMP_DIR='tempdir'
OUT_DIR='output'

k=24
nk1=$k
nk2=$k
nk3=$k

# species.in -------------------------------------------------------------------
cat > species.in << EOF
ATOMIC_SPECIES
Nb  92.906  Nb.upf

ATOMIC_POSITIONS crystal
Nb       0.000000000   0.000000000   0.000000000

EOF


lfail=false

OCP="-rv"
ORM="-vrf"
OMV="-v"
OTAR="-Jcvf"
ORSYNC="-avz --delete"

file_check() {
  if [[ -f "$1" ]]; then
    $ECHO "$1 exists."
  else
    $ECHO "Error: $1 does not exist."
    lfail=true
  fi
}

dir_check() {
  if [[ -d "$1" ]]; then
    $ECHO "$1 exists."
  else
    $ECHO "Error: $1 does not exist."
    lfail=true
  fi
}

ngauss=1
if [[ "$SMEARING" == "mp" ]]
then
  ngauss=1
elif [[ "$SMEARING" == "gauss" ]]
then
  ngauss=0
elif [[ "$SMEARING" == "fd" ]]
then
  ngauss="-99"
fi

$ECHO "smearing: " $SMEARING " k: " $k




#. clean.sh
mkdir -p $OUT_DIR
mkdir -p $TMP_DIR

# PREPARE ======================================================================

# check input ------------------------------------------------------------------
$ECHO "lrun_scf = " $lrun_scf
$ECHO "lrun_rmta = " $lrun_rmta
$ECHO "ncpu = " $ncpu
$ECHO "occupations = " $OCCUPATIONS
$ECHO "smearing = " $SMEARING
$ECHO "degauss = " $DEGAUSS
$ECHO "ngauss = " $ngauss
$ECHO "k = " ${k}
$ECHO "IBRAV = " $IBRAV
$ECHO "ALAT = " $ALAT
$ECHO "ECUT = " $ECUT
$ECHO "NAT = " $NAT
$ECHO "NTYP = " $NTYP
$ECHO ""

# alat.in ----------------------------------------------------------------------
echo $ALAT > alat.in



# SCF ==========================================================================
if [[ "$lrun_scf" == "true" && "$lfail" != "true" ]]
then
TASK='scf'
SUFFIX=$TASK

$ECHO ""
$ECHO "$SUFFIX"
$ECHO "-------------->"

NAME=$PREFIX.$TASK
cat > $NAME.in << EOF
&control
 calculation = '$TASK',
 restart_mode = 'from_scratch',
 prefix = '$PREFIX',
 tprnfor = .true.,
 tstress =.true.,
 pseudo_dir = '$PSEUDO_DIR',
 outdir = '$TMP_DIR'
 etot_conv_thr = 1.0d-5
 forc_conv_thr = 1.0d-4
/

&system
 ibrav = $IBRAV
 celldm(1) = $ALAT
 nat = $NAT
 ntyp = $NTYP
 ecutwfc = $ECUT
 occupations = '$OCCUPATIONS'
 degauss = $DEGAUSS
 smearing = '$SMEARING'
/

&electrons
 diagonalization = 'david'
 mixing_mode = 'plain'
 mixing_beta = 0.7
 conv_thr =  1.0d-12
/

&cell
  press = $PRESS ! kbar
/

K_POINTS AUTOMATIC
$nk1 $nk2 $nk3 0 0 0

EOF

name_check="species.in"
file_check $name_check
cat $name_check >> $NAME.in

$ECHO "  running the scf calculation for $PREFIX..."
$EXEC -n $ncpu $BIN_DIR/pw.x < $NAME.in > $OUT_DIR/$NAME.out
$ECHO "$SUFFIX is done"

if [ $? -ne 0 ]; then
  $ECHO “Error: Failed $SUFFIX.”
  lfail=true
fi

mkdir -p tempdir_$SUFFIX
rsync $ORSYNC tempdir tempdir_$SUFFIX/

if [[ $lsave_tmp_dir_tar == "true" ]]
then
  XZ_OPT="-T${ncpu}" tar $OTAR ${SUFFIX}_out.tar.xz $TMP_DIR
fi

grep "Fermi" output/$NAME.out \
| sed -e 's/the Fermi energy is//g' \
| sed -e 's/ev//g' > efermi.in


fi





# RMTA =========================================================================
if [[ "$lrun_rmta" == "true" && "$lfail" != "true" ]]
then
TASK='rmta'
SUFFIX=$TASK

$ECHO ""
$ECHO "$SUFFIX"
$ECHO "-------------->"

name_check="${TMP_DIR}_scf/${TMP_DIR}"
if [[ -d $name_check ]]
then
dir_check $name_check
rsync $ORSYNC $name_check .
else
name_check=scf_out.tar.xz
file_check $name_check
rm $ORM ${TMP_DIR}
tar -xvf $name_check
mkdir -p ${TMP_DIR}_scf
rsync $ORSYNC ${TMP_DIR} ${TMP_DIR}_scf/
fi


NAME=$PREFIX.$SUFFIX
$ECHO "  running RMTA for $PREFIX..."
cat > ${NAME}.in << EOF
&rmta
  prefix = '$PREFIX'
  outdir = '$TMP_DIR'
  lwrite_dat = .true.
/
EOF
$EXEC -n 1 $BIN_DIR/rmta.x < $NAME.in > $OUT_DIR/$NAME.out

$ECHO "$SUFFIX is done"

if [ $? -ne 0 ]
then
  $ECHO “Error: Failed $SUFFIX.”
  lfail=true
fi

# rm $ORM $TMP_DIR

fi




# SUMMARY ======================================================================
if [[ "$lfail" == "true" ]]
then
$ECHO "FAIL"
else

if [[ $lsave_tmp_dir != "true" ]]
then
  rm $ORM ${TMP_DIR}
  rm $ORM ${TMP_DIR}_scf
fi

$ECHO "DONE"
fi
