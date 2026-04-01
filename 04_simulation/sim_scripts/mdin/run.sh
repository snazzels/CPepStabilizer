set -e

rootdir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDP=$rootdir/mdin

# rst + top (relative to rootdir)
top=system_wb.prmtop
start=system_wb.rst

#-----------------------------------------------------------------------
#EM
mkdir -p $rootdir/em/
pmemd.cuda -O -i $MDP/em_1.in \
              -p $rootdir/${top} \
              -c $rootdir/${start} \
              -ref $rootdir/${start} \
	      -o $rootdir/em/em_1.out \
              -r $rootdir/em/em_1.rst \
              -inf $rootdir/em/em_1.mdinfo

pmemd.cuda -O -i $MDP/em_2.in \
              -p $rootdir/${top} \
              -c $rootdir/em/em_1.rst \
	      -ref $rootdir/em/em_1.rst \
              -o $rootdir/em/em_2.out \
              -r $rootdir/em/em_2.rst \
              -inf $rootdir/em/em_2.mdinfo

echo "EM done"


#-----------------------------------------------------------------------
#HEATING (CONTINUOUS)
mkdir -p $rootdir/heat/
pmemd.cuda -O -i "$MDP"/heat.in \
              -p $rootdir/${top}\
              -c $rootdir/em/em_2.rst \
              -ref $rootdir/em/em_2.rst \
              -o $rootdir/heat/heat.out \
              -r $rootdir/heat/heat.rst \
              -x $rootdir/heat/heat.nc \
              -inf $rootdir/heat/heat.mdinfo

echo "Heating done"

#-----------------------------------------------------------------------
#EQUILIBRATION
#NVT
mkdir -p $rootdir/nvt/
pmemd.cuda -O -i "$MDP"/nvt_1.in \
              -p $rootdir/${top} \
              -c $rootdir/heat/heat.rst \
              -ref $rootdir/heat/heat.rst \
              -o $rootdir/nvt/nvt_1.out \
              -r $rootdir/nvt/nvt_1.rst \
              -x $rootdir/nvt/nvt_1.nc \
              -inf $rootdir/nvt/nvt_1.mdinfo

pmemd.cuda -O -i "$MDP"/nvt_2.in \
              -p $rootdir/${top} \
              -c $rootdir/nvt/nvt_1.rst \
              -ref $rootdir/nvt/nvt_1.rst \
              -o $rootdir/nvt/nvt_2.out \
              -r $rootdir/nvt/nvt_2.rst \
              -x $rootdir/nvt/nvt_2.nc \
              -inf $rootdir/nvt/nvt_2.mdinfo

echo "NVT equilibration done"


#NPT
mkdir -p $rootdir/npt/
pmemd.cuda -O -i "$MDP"/npt_1.in \
              -p $rootdir/${top} \
              -c $rootdir/nvt/nvt_2.rst \
              -ref $rootdir/nvt/nvt_2.rst \
              -o $rootdir/npt/npt_1.out \
              -r $rootdir/npt/npt_1.rst \
              -x $rootdir/npt/npt_1.nc \
              -inf $rootdir/npt/npt_1.mdinfo

pmemd.cuda -O -i "$MDP"/npt_2.in \
              -p $rootdir/${top} \
              -c $rootdir/npt/npt_1.rst \
              -ref $rootdir/npt/npt_1.rst \
              -o $rootdir/npt/npt_2.out \
              -r $rootdir/npt/npt_2.rst \
              -x $rootdir/npt/npt_2.nc \
              -inf $rootdir/npt/npt_2.mdinfo

pmemd.cuda -O -i "$MDP"/npt_3.in \
              -p $rootdir/${top} \
              -c $rootdir/npt/npt_2.rst \
              -ref $rootdir/npt/npt_2.rst \
              -o $rootdir/npt/npt_3.out \
              -r $rootdir/npt/npt_3.rst \
              -x $rootdir/npt/npt_3.nc \
              -inf $rootdir/npt/npt_3.mdinfo

echo "NPT equilibration done"

mkdir -p $rootdir/free_run/
pmemd.cuda -O -i "$MDP"/run_free.in \
              -p $rootdir/${top} \
              -c $rootdir/npt/npt_3.rst \
              -o $rootdir/free_run/run.out \
              -r $rootdir/free_run/run.rst \
              -x $rootdir/free_run/run.nc \
              -inf $rootdir/free_run/run.mdinfo
