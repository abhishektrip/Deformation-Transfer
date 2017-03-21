#include <iostream>
#include <string>
#include "MeshConstraints.h"

using namespace DeformationTransfer::Correspondence;
using namespace std;

MeshConstraints::MeshConstraints()
{
}


MeshConstraints::~MeshConstraints()
{
}

bool MeshConstraints::ReadConstraintFile(string path)
{
	FILE* file = fopen(path.c_str(), "r");
	char prefix[3];   // prefix reader - 3 chars should be sufficient. 
	bool success = false;
	if (file == nullptr)
	{
		cout << "Could not open file :" << path << endl;		
	}
	else
	{		
		fscanf(file, "%d", &numConstraints);
		if (numConstraints > 0)
		{
			constrainedList.reserve(numConstraints);

			for (int i = 0; i < numConstraints; i++)
			{
				int key, value;
				fscanf(file, "%d,%d", &key, &value);
				constrainedList.insert({ key,value });
			}
			success = true;
		}		
	}
	// Close file and return status. 
	fclose(file);
	return success;
}
