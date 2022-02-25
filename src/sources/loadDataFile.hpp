#ifndef LOADDATAFILE_HPP
#define LOADDATAFILE_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>
#include <set>

using namespace std;

template <typename T>
void loadDataFile(const string& fileName, const set<T>& sections, map<T, vector<string> >& dataFile) {
  ifstream f(fileName.c_str());
  if (! f.is_open()){
    cout << "Unable to open file " << fileName << " ." << endl;
    exit(1);
  }

  string loadPart;
  T partTest;
  T oneSection;
  typename set<T>::const_iterator it;

  while (!f.eof()) {
    stringstream part;
    f >> loadPart;
    part << loadPart;
    part >> partTest;
    it = sections.find(partTest);
    
    string comment;
    comment.append(loadPart, 0, 2);
    
    if (comment == "//") {
      string radek;
      getline(f, radek);
    }
    else if (comment == "/*") {
      string endComment;
      endComment.append(loadPart, loadPart.size()-2, 2);
     
      if (endComment != "*/") {
	string trash;
	while(!f.eof() && endComment != "*/") {
	  endComment.clear();
	  f >> trash;
	  if (trash.size() > 1) {
	    endComment.append(trash, trash.size()-2, 2);
	  }
	  else {
	    endComment = trash;
	  }
	}
      }
    }
    else if (it != sections.end()) {
      oneSection = partTest;
    }
    else {
      dataFile[oneSection].push_back(loadPart);
    }
  }

  f.close();
};


#endif
