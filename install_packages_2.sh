wget https://github.com/jgm/pandoc/releases/download/2.9.2.1/pandoc-2.9.2.1-1-amd64.deb
sudo dpkg -i pandoc-2.9.2.1-1-amd64.deb
apt-get install -y libxml2-dev
apt-get install -y libcurl4-openssl-dev
Rscript install_packages_2.R
