#ifndef FIELD2_HPP
#define FIELD2_HPP

#include <iostream>
#include <cstdlib>

using namespace std;

template <typename T>
class Field2 {
  int imin;
  int imax;
  int jmin;
  int jmax;

  T *data;

  bool allocated;

public:
  Field2(): imin(0), imax(0), jmin(0), jmax(0) {allocated = false;};

  Field2(const Field2& fl2) {
    allocated = false;
    allocate(fl2.Imin(), fl2.Imax(), fl2.Jmin(), fl2.Jmax());

    int size_i = imax - imin;
    int size_j = jmax - jmin;
    int totalSize = size_i * size_j;
    
    for (int i=0; i<totalSize; i++) {
      data[i] = fl2.data[i];
    }
  };

  Field2(int Imin, int Imax, int Jmin, int Jmax): imin(Imin), imax(Imax), jmin(Jmin), jmax(Jmax) {
    //alokovani pameti
    int size_i = imax - imin;
    int size_j = jmax - jmin;
    int totalSize = size_i * size_j;

    data = new T[totalSize];
    

    allocated = true;
  };

  ~Field2() {
    if (allocated) {
      delete [] data;
    }
  };

  void allocate(const int Imin, const int Imax, const int Jmin, const int Jmax) {
    if (allocated) {
      cout << "Field2 uz je naalokovane, pouzijte fci reallocate()!" << endl;
      exit(0);
    }

    imin = Imin;
    imax = Imax;
    jmin = Jmin;
    jmax = Jmax;

    //alokovani pameti
    int size_i = imax - imin;   //pocet prvku ve smeru i
    int size_j = jmax - jmin;   //pocet prvku ve smeru j
    int totalSize = size_i * size_j;

    data = new T[totalSize];

    allocated = true;
  }

  void reallocate(const int Imin, const int Imax, const int Jmin, const int Jmax) {
    if (!allocated) {
      cout << "Field2 jeste neni naalokovano, pouzijte fci allocate()!" << endl;
      exit(1);
    }

    //dealokovani pameti
    delete [] data;

    allocated = false;

    allocate(Imin, Imax, Jmin, Jmax);
  }

  int Imin() const {return imin;}
  int Imax() const {return imax;}
  int Jmin() const {return jmin;}
  int Jmax() const {return jmax;}

  inline T* operator[](const int i) const {
    int size_j = jmax - jmin;
    return data + (i-imin)*size_j - jmin;
  }

  Field2& operator=(const Field2& fl2) {
    if (allocated) {
      reallocate(fl2.Imin(), fl2.Imax(), fl2.Jmin(), fl2.Jmax());
    }
    else {
      allocate(fl2.Imin(), fl2.Imax(), fl2.Jmin(), fl2.Jmax());
    }

    int size_i = imax - imin;
    int size_j = jmax - jmin;
    int totalSize = size_i * size_j;
    
    for (int i=0; i<totalSize; i++) {
      data[i] = fl2.data[i];
    }

    return *this;
  }
};

#endif
