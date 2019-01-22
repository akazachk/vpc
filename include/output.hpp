#pragma once
#include <string>
#include <cstdio>

/***********************************************************************/
/**
 * @brief Writes error and closes file myfile.
 */
inline void writeErrorToLog(std::string text, FILE *myfile) {
	if (myfile == NULL || myfile == stdout)
		return;
	fprintf(myfile, "||%c%s", ',', text.c_str());
	fclose(myfile);
}
