// veb_Dominik_Dragon.cpp: Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cstdlib>
#include <cmath>
#include "veb.cpp"

void test1() {
	int64_t res;
	//int64_t universe = 16777216; // 2 ^ 24 (16 млн)
	int64_t universe = 4294967296; // 2 ^32 (4294 млн)
	TvEB * tree = new TvEB(universe);

	res = vEB_insert(tree, 2);

	if (!res) printf("insert problem veb 2_1\n");
	// Дуликаты не поддерживаются
	res = vEB_insert(tree, 2);
	if (!res) printf("insert problem veb 2_2\n");

	res = vEB_insert(tree, 249236743);
	if (!res) printf("insert problem veb 249236743\n");

	res = vEB_insert(tree, 49236743);
	if (!res) printf("insert problem veb 49236743\n");

	res = vEB_insert(tree, 236743);
	if (!res) printf("insert problem veb 236743\n");

	res = vEB_insert(tree, 743);
	if (!res) printf("insert problem veb 743\n");

	int64_t succ;
	res = vEB_succ(tree, 3, succ);
	int64_t pred;
	res = vEB_pred(tree, 13, pred);

	res = vEB_find(tree, 8);

	res = vEB_delete(tree, 8);

	if (vEB_min(tree, res)) {
		//printf("minimum =%d\n", res);
		std::cout << "minimum =" << res << std::endl;
	}
	if (vEB_max(tree, res)) {
		//printf("maximum =%d\n", res);
		std::cout << "maximum =" << res << std::endl;
	}
}

int64_t N=64;
int64_t* a = new int64_t[N];

int main()
{
	int64_t res;
	int64_t universe = 4294967296; // 2 ^32 (4294 млн)
	TvEB * tree = new TvEB(universe);

	//test1();
	for (int64_t i = 0; i < N; i++) {
		a[i] = abs(rand())+1;
		res = vEB_insert(tree, a[i]);
	}

	int64_t ic21 = 1;
	if (vEB_min(tree, res)) printf("%lld %lld \n", ic21, res);
	ic21++;
	int64_t succ;
	while (vEB_succ(tree, res, succ)) {
		res = succ;
		printf("%lld %lld \n",ic21, res);
		ic21++;
	}
	getchar();
	
    return 0;
}

