# Dependencies #

The project uses cmake to build the program. The following command can be used to install cmake:

  * sudo apt-get install cmake

The usage of cmake allows the project to enable the opencv Qt qui. Qt5 can be installed using the follwing commands:


  * sudo apt-add-repository ppa:ubuntu-sdk-team/ppa
  * sudo apt-get update
  * sudo apt-get install qtdeclarative5-dev/

Finally the image processing library used in this project is opencv, which can be installed using the following instructions:

  * http://milq.github.io/install-opencv-ubuntu-debian/

After all those dependencies are installed, the project can be built by going to the project's directory and running the following commands:

  * cmake CMakeList.txt
  * make

Finally the program can be executed by running with the following command:

  * /Method path\_to\_image