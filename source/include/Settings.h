#ifndef _Settings_h_
#define _Settings_h_

class Settings{

 public:

  Settings();
  ~Settings();

  // Getters
  void getR0();

  // Setters
  void setR0(double r0){R0=r0;}

  // Print the settings
  void printSettings();


 private:

  double R0;

};


#endif
