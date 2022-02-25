#ifndef FINDSECTION_HPP
#define FINDSECTION_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

template <typename T, typename S>
void findSection(const map<T, vector<string> >& dataFile, const string& name, const T& section, S& value) {
  string dataName;
  stringstream streamValue;

  typename map<T, vector<string> >::const_iterator it;
  it = dataFile.find(section);

  if (it != dataFile.end()) {
    const vector<string>& dataSection = it->second;
    for (int i=0; i<dataSection.size(); i++) {
      dataName = dataSection[i];
      if (dataName == name) {
	streamValue << dataSection[i+1];
	streamValue >> value;
	goto end;
      }
    }
    cout << "Not found " << name << " in a section " << section << "!" << endl;
    exit(1);
  }
  else {
    cout << "Not found section " << section << "!" << endl;
    exit(1);
  }

 end:;
}

#endif
