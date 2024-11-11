# openmm-mdrun
A wrapper on openmm for basic mdrun

## 1. Installation
### 1.1 Environment
Please use `conda` or `mamba` to create a environment.  
```bash
mamba create -n open-env -c conda-forge openmm openmmtools mdtraj numpy pandas scipy mdanalysis pymbar-core tqdm
```

### 1.2 Install
```bash
mamba activate open-env
git clone https://github.com/huichenggong/openmm-mdrun.git
cd openmm-mdrun
pip install .
```

## 2. Usage
Several command line tools are provided.  

### 2.1 `openmm_mdrun`
This tool is for basic mdrun.  
```bash
openmm_mdrun \
        -sys  system.xml.gz \
        -p    step5_input.psf \
        -mdp  npt.mdp \
        -t    ../01-nvt/npt.xml \
        -cpt  npt.xml \
        -deffnm npt \
        -maxh 12
```

### 2.2 `openmm_GCNCMC`
This tool is for GCNCMC (Grand Canonical Non-Equilibrium Candidate Monte Carlo).  
```bash
i_old=0
for i in {1..10} ; do
    if [ ! -f md${i}.rst7 ] ; then
        echo "Starting md$i"
        openmm_GCNCMC \
            -sys ../../../../01-EM/system.xml.gz \
            -p   ../../../../step5_input.psf \
            -mdp ../npt.mdp \
            -t ../../02-noRes/state_119_ghosts_water30.rst7 \
            -rst     md${i_old}.rst7 \
            -ighosts md${i_old}.dat \
            -ilog    md${i_old}.log \
            -atom    ../atom.json \
            -deffnm md$i \
            -maxh 12 --debug
        break
    else
        echo "md${i}.rst7 exists"
    fi
    i_old=$i
done
```