#include "canny.h"
using namespace std;

int main() {
	canny bigben("./test/lena.jpg");
	bigben.run();
	return 0;
}