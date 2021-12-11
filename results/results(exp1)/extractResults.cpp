#include <iostream>
#include <fstream>
#include <dirent.h>
#include <vector>
#include <map>
#include <string>
#include <sstream>

using namespace std;

vector<string> getFilesList() {
	vector<string> filesList;

	DIR *dr;
	struct dirent *en;
	dr = opendir("."); //open all or present directory
	if (dr) {
		while ((en = readdir(dr)) != NULL) {
			//printf("%s\n", en->d_name); //print all directory name
			filesList.push_back(en->d_name);
		}
		closedir(dr); //close all directory
	}
	return filesList;
}

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

int main() {
	std::ofstream results;
	std::ifstream inputFile;

	map<string, vector<string>> table;

	vector<string> files = getFilesList();
	for(string s : files) {
		cout << s << endl;
		
		if (!s.compare("."))
			continue;

		else if (!s.compare(".."))
			continue;

		else if (!s.compare("extractResults.cpp"))
			continue;

		else if (!s.compare("extract_results"))
			continue;

		else if (!s.compare("comparsion.csv"))
			continue;

		else if (!s.compare("comparsion_2.csv"))
			continue;

		table["config"].push_back(s);
		inputFile.open(s);

		string line;
		while (getline(inputFile, line)) {
			vector<string> words = split(line, ';');
			table[words[0]].push_back(words[3]);
		}

		inputFile.close();
	}


	results.open("comparsion.csv");
	for (auto const &row : table) {
		results << row.first;
		for (string s : row.second) {
			results << ";" << s;
		}
		results << endl;
	}
	
	results.close();
	return 0;
}