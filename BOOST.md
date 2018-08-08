# How to install Boost unit test framework library

We are using the boost library for our c++ testing framework following the 1.62 and onward syntax. There are several solutions to install it.


## ubuntu/debian system wide
```
sudo apt-get install libboost-test-dev
```
Depending on your version of these distributions, you might not get a recent enough verison of the boost library.

## Mac
It seems that the brew recepy of boost is broken (compilation errors). Custom or Conda installations are advised. 

## Custom installation
We use 7z to unpack (sudo apt-get install p7zip-full to install on ubuntu)
```
export BOOST_DIR=/path/to/download/directory
wget https://sourceforge.net/projects/boost/files/boost/1.67.0/boost_1_67_0.7z/download -O $BOOST_DIR/boost_1_67_0.7z
7z x $BOOST_DIR/boost_1_67_0.7z
cd $BOOST_DIR/boost_1_67_0
mkdir boost_root
./bootstrap.sh --with-libraries=test --prefix=$BOOST_DIR/boost_root
./b2 install
export Boost_NO_SYSTEM_PATHS=ON
export BOOST_ROOT=$BOOST_DIR/boost_root
```
Cmake FindBoost might not find this installation of Boost Test so we set Boost_NO_SYSTEM_PATHS to avoid finding another installation of Boost and we set BOOST_ROOT to this new installation.

## Using Conda
Conda distributes Boost library in the conda forge channel.
+ First install miniconda if not already available:
```
export MINICONDA_DOWNLOAD_DIR=/path/to/download/directory
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $MINICONDA_DOWNLOAD_DIR/miniconda.sh
bash $MINICONDA_DOWNLOAD_DIR/miniconda.sh
bash $MINICONDA_DOWNLOAD_DIR/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH" 
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc
source $HOME/.bashrc
conda create -n myenv python=3.6
```

+ Install boost from conda forge
```
source activate myenv
conda install -c conda-forge boost
export Boost_NO_SYSTEM_PATHS=ON 
export BOOST_ROOT=/path/to/myenv
```
Following the previous example, /path/to/myenv -> $HOME/miniconda/envs/myenv.

Cmake FindBoost might not find this installation of Boost so we set Boost_NO_SYSTEM_PATHS to avoid finding another installation of Boost and we set BOOST_ROOT to this new installation.