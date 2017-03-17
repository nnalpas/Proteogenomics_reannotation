


# Create software location directory
mkdir ~/Software/blast+

# Go to software repository
cd Software/

# Download the software
wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-src.zip'
  
# Load the java module
module load devel/java_jdk/1.8.0u112

# Unzip installation files using java (unzip command was not installed at the time)
jar -xvf ncbi-blast-2.6.0+-src.zip

# Make executable all files in software installation directory
chmod 770 -R ncbi-blast-2.6.0+-src

# Move to installation directory
cd ncbi-blast-2.6.0+-src/c++/

# Configure the installation
./configure

# Move to build directory
cd ReleaseMT/build/

# Install software
make all_r

# Rename software directory
mv ~/Software/ncbi-blast-2.6.0+-src/ ~/Software/blast+/2.6.0


