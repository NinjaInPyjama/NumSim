FILE* handle = fopen(file,"r");
double inval[2];
// Code-Schnipsel zum parsen von Parameter/Geometrie Files

#include <cstdio>
#include <cstring>
#include <cstdlib>

//--------------------------------------------------------------------------------------------------

char name[20];
while (!feof(handle)) {
	if (!fscanf(handle, "%s =", name)) continue;
	if (strcmp(name,"size") == 0) {
		if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
			_size[0] = inval[0];
			_size[1] = inval[1];
		}
		continue;
	}
	if (strcmp(name,"length") == 0) {
		if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
			_length[0] = inval[0];
			_length[1] = inval[1];
		}
		continue;
	}
	if (strcmp(name,"velocity") == 0) {
		if (fscanf(handle," %lf %lf\n",&inval[0],&inval[1])) {
			_velocity[0] = inval[0];
			_velocity[1] = inval[1];
		}
		continue;
	}
	if (strcmp(name,"pressure") == 0) {
		if (fscanf(handle," %lf\n",&inval[0]))
			_pressure = inval[0];
		continue;
	}
}
fclose(handle);

//--------------------------------------------------------------------------------------------------

FILE* handle = fopen(file,"r");
double inval;
char name[20];
while (!feof(handle)) {
	if (!fscanf(handle, "%s = %lf\n", name, &inval)) continue;
	if (strcmp(name,"re") == 0) _re = inval;
	else if (strcmp(name,"omg") == 0) _omega = inval;
	else if (strcmp(name,"alpha") == 0) _alpha = inval;
	else if (strcmp(name,"dt") == 0) _dt = inval;
	else if (strcmp(name,"tend") == 0) _tend = inval;
	else if (strcmp(name,"iter") == 0) _itermax = inval;
	else if (strcmp(name,"eps") == 0) _eps = inval;
	else if (strcmp(name,"tau") == 0) _tau = inval;
	else printf("Unknown parameter %s\n",name);
}
fclose(handle);
