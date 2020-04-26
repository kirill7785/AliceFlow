/*******************************************************************************
* Copyright (C) 2016 Dominik Dragoun                                           *
*                                                                              *
* This file is part of my Van Emde Boas tree data structure implementation in  *
* C/C++. It is released under MIT License, which should be distributed with    *
* this file. It is also avaible at <https://opensource.org/licenses/MIT>.      *
*******************************************************************************/

/***************************************************************************//**
																			 * @file veb.cpp
																			 *
																			 * @brief      File containing definition of a class implementing Van Emde Boas
																			 *             tree data structure.
																			 * @author     Dominik Dragoun (dominik@dragoun.com)
																			 * @date       June, 2016
																			 * @copyright  Copyright (C) 2016 Dominik Dragoun.
																			 * @license    This project is released undes the MIT License.
																			 ******************************************************************************/

#pragma once
#ifndef DOMINIK_DRAGON_VEB_CPP
#define DOMINIK_DRAGON_VEB_CPP 1

#include "veb.h"

// Исправлено 21.03.2019
// Был конфликт имён члена класса min с функцией min(a,b) языка СИ.
// Аналогично для поля класса max. Поля классов min и max переименованы 
// в данные класса veb_min, veb_max. Теперь конфликт имён отсутсвует.
TvEB::TvEB(int64_t uniSize)
	: uni(powTwoRoundUp(uniSize)), uniSqrt((int64_t)(sqrt(uni))),
	lowerUniSqrt((int64_t)(lowerSqrt(uni))), higherUniSqrt((int64_t)(higherSqrt(uni))),
	veb_min(UNDEFINED), veb_max(UNDEFINED), summary(NULL)
{
	if (uniSize <= 0)
	{
		std::cerr << "universe size of TvEB must be bigger than 0" << std::endl;
		return;
	}

	if (uni > 2)
	{
		cluster = new TvEB *[higherUniSqrt];
		for (int64_t i = 0; i < higherUniSqrt; ++i)
		{
			cluster[i] = NULL;
		}
	}
	else
	{
		cluster = NULL;
	}
}

TvEB::~TvEB()
{
	if (summary) delete summary;
	summary = NULL;
	if (cluster)
	{
		for (int64_t i = 0; i < higherUniSqrt; ++i)
		{
			if (cluster[i]) delete cluster[i];
		}
		delete[] cluster;
	}
	cluster = NULL;
}

int64_t powTwoRoundUp(int64_t x)
{
	if (x < 1) return 0;
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
}

double lowerSqrt(int64_t val)
{
	return pow(2, floor(log2(val) / 2));
}

double higherSqrt(int64_t val)
{
	return pow(2, ceil(log2(val) / 2));
}

int64_t low(TvEB * tree, int64_t val)
{
	return val % (int64_t)lowerSqrt(tree->uni);
}

int64_t high(TvEB * tree, int64_t val)
{
	return (int64_t)(val / lowerSqrt(tree->uni));
}

int64_t index(TvEB * tree, int64_t high, int64_t low)
{
	return (int64_t)(high * lowerSqrt(tree->uni)) + low;
}

bool vEB_min(TvEB * tree, int64_t & res)
{
	if (tree)
	{
		res = tree->veb_min;
		return true;
	}
	return false;
}

bool vEB_max(TvEB * tree, int64_t & res)
{
	if (tree)
	{
		res = tree->veb_max;
		return true;
	}
	return false;
}

bool vEB_insert(TvEB *& tree, int64_t val, int64_t parentUniSqrt)
{
	if (!tree)
	{
		tree = new TvEB(parentUniSqrt);
	}

#ifdef DEBUG
	DEBUG_OS << "inserting " << val << " to tree " << tree
		<< " of size " << tree->uni << DEBUG_OS_ENDL;
#endif /* DEBUG */

	if (val < 0 || val >= tree->uni) return false;

	if (tree->veb_min == val || tree->veb_max == val) return false;

	if (tree->veb_min == UNDEFINED)
	{
		tree->veb_min = tree->veb_max = val;
		return true;
	}

	if (val < tree->veb_min)
	{
		int64_t tmp = val;
		val = tree->veb_min;
		tree->veb_min = tmp;
	}

	if (val > tree->veb_max)
	{
		tree->veb_max = val;
	}

	if (tree->uni > 2)
	{
		int64_t lowVal = low(tree, val);
		int64_t highVal = high(tree, val);
		if (!tree->cluster[highVal])
		{
			if (!vEB_insert(tree->summary, highVal, tree->higherUniSqrt)) return false;
		}

		if (!vEB_insert(tree->cluster[highVal], lowVal, tree->lowerUniSqrt)) return false;
	}
	return true;
}

bool vEB_delete(TvEB *& tree, int64_t val)
{
	if (!tree) return false;

#ifdef DEBUG
	DEBUG_OS << "deleting " << val << " from tree " << tree
		<< " of size " << tree->uni << DEBUG_OS_ENDL;
#endif /* DEBUG */

	if (val < 0 || val >= tree->uni) return false;
	if (tree->veb_min > val || tree->veb_max < val) return false;

	if (tree->veb_min == val)
	{
		int64_t i;
		if (!vEB_min(tree->summary, i) || i == UNDEFINED)
		{
			if (tree->veb_min != tree->veb_max)
			{
				tree->veb_min = tree->veb_max;
				return true;
			}

			tree->veb_min = tree->veb_max = UNDEFINED;
			delete tree;
			tree = NULL;
			return true;
		}

		val = tree->veb_min = index(tree, i, tree->cluster[i]->veb_min);
	}

	if (tree->uni > 2)
	{
		int64_t highVal = high(tree, val);
		if (!vEB_delete(tree->cluster[highVal], low(tree, val))) return false;

		int64_t tmp;
		if (!vEB_min(tree->cluster[highVal], tmp) || tmp == UNDEFINED)
		{
			if (!vEB_delete(tree->summary, highVal)) return false;
		}
	}

	if (tree->veb_max == val)
	{
		int64_t tmp;
		if (!vEB_max(tree->summary, tmp) || tmp == UNDEFINED)
		{
			tree->veb_max = tree->veb_min;
		}
		else
		{
			int64_t i;
			if (!vEB_max(tree->summary, i)) return false;
			tree->veb_max = index(tree, i, tree->cluster[i]->veb_max);
		}
	}
	return true;
}

bool vEB_find(TvEB * tree, int64_t val)
{
	if (!tree) return false;

#ifdef DEBUG
	DEBUG_OS << "looking for " << val << " in tree " << tree << DEBUG_OS_ENDL;
#endif /* DEBUG */

	if (val < 0 || val >= tree->uni) return false;
	if (tree->veb_min > val || tree->veb_max < val) return false;
	if (tree->veb_min == val) return true;
	if (!tree->summary)
	{
		return tree->veb_max == val;
	}
	if (!vEB_find(tree->cluster[high(tree, val)], low(tree, val)))
		return false;
	return true;
}

bool vEB_succ(TvEB * tree, int64_t val, int64_t & res)
{
	if (!tree) return false;

#ifdef DEBUG
	DEBUG_OS << "looking for successor of " << val << " in tree " << tree
		<< " of size " << tree->uni << DEBUG_OS_ENDL;
#endif /* DEBUG */

	if (val < -1 || val >= tree->uni) return false;

	if (tree->veb_min > val)
	{
		res = tree->veb_min;
		return true;
	}

	if (!tree->summary)
	{
		if (tree->veb_max > val)
		{
			res = tree->veb_max;
			return true;
		}
		return false;
	}

	int64_t lowVal = low(tree, val);
	int64_t highVal = high(tree, val);
	int64_t i = highVal;
	int64_t j = UNDEFINED;
	int64_t tmp;
	if (vEB_max(tree->cluster[i], tmp) && lowVal < tmp)
	{
		if (!vEB_succ(tree->cluster[i], lowVal, j)) return false;
	}
	else
	{
		if (!vEB_succ(tree->summary, highVal, i))
		{
			if (tree->veb_max > val)
			{
				res = tree->veb_max;
				return true;
			}
			return false;
		}
		if (!vEB_min(tree->cluster[i], j)) return false;
	}

	res = index(tree, i, j);
	return true;
}

bool vEB_pred(TvEB * tree, int64_t val, int64_t & res)
{
	if (!tree) return false;

#ifdef DEBUG
	DEBUG_OS << "looking for predecessor of " << val << " in tree " << tree
		<< " of size " << tree->uni << DEBUG_OS_ENDL;
#endif /* DEBUG */

	if (val < 0 || val > tree->uni) return false;

	if (tree->veb_max < val)
	{
		res = tree->veb_max;
		return true;
	}

	if (!tree->summary)
	{
		if (tree->veb_min < val)
		{
			res = tree->veb_min;
			return true;
		}
		return false;
	}

	int64_t lowVal = low(tree, val);
	int64_t highVal = high(tree, val);
	int64_t i = highVal;
	int64_t j = UNDEFINED;
	int64_t tmp;
	if (vEB_min(tree->cluster[i], tmp) && lowVal > tmp)
	{
		if (!vEB_pred(tree->cluster[i], lowVal, j)) return false;
	}
	else
	{
		if (!vEB_pred(tree->summary, highVal, i))
		{
			if (tree->veb_min < val)
			{
				res = tree->veb_min;
				return true;
			}
			return false;
		}
		if (!vEB_max(tree->cluster[i], j)) return false;
	}

	res = index(tree, i, j);
	return true;
}

void vEB_print(TvEB * tree, std::ostream & os)
{
	if (!tree) return;
	os << "tree: " << tree << std::endl;
	os << "min: " << tree->veb_min << ", max: " << tree->veb_max << std::endl;
	os << "uni: " << tree->uni << ", uniSqrt: " << tree->uniSqrt << std::endl;
	os << "lowerUniSqrt: " << tree->lowerUniSqrt;
	os << ", higherUniSqrt: " << tree->higherUniSqrt << std::endl;
	os << "summary: " << tree->summary << std::endl;
	if (tree->uni > 2)
	{
		for (int64_t i = 0; i < tree->higherUniSqrt; ++i)
		{
			os << "cluster " << i << ": " << tree->cluster[i] << std::endl;
		}
	}
	else
	{
		os << "cluster " << tree->cluster << std::endl;
	}
}

#endif /*DOMINIK_DRAGON_VEB_CPP*/