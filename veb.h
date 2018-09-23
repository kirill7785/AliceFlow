
/*******************************************************************************
* Copyright (C) 2016 Dominik Dragoun                                           *
*                                                                              *
* This file is part of my Van Emde Boas tree data structure implementation in  *
* C/C++. It is released under MIT License, which should be distributed with    *
* this file. It is also avaible at <https://opensource.org/licenses/MIT>.      *
*******************************************************************************/

/***************************************************************************//**
																			 * @file veb.h
																			 *
																			 * @brief      File containing declarations of a class implementing Van Emde
																			 *             Boas tree data structure.
																			 * @author     Dominik Dragoun (dominik@dragoun.com)
																			 * @date       June, 2016
																			 * @copyright  Copyright (C) 2016 Dominik Dragoun.
																			 * @license    This project is released undes the MIT License.
																			 ******************************************************************************/

#ifndef __VEB_H_458976543568798867538649758687654752463747856374562543646__
#define __VEB_H_458976543568798867538649758687654752463747856374562543646__

#include <cmath>
#include <climits>
#include <cstdlib>
#include <cstdio>
#include <iostream>

																			 // #define DEBUG
#define DEBUG_OS std::cout
#define DEBUG_OS_ENDL std::endl

#define UNDEFINED INT_MIN

																			 /***************************************************************************//**
																																						  * @brief      Struct containing the Van Emde Boas tree.
																																						  *
																																						  * @details    It is a data structure which implements an associative array with
																																						  *             m-bit integer keys. It supports operations of ordered associative
																																						  *             array (Insert, Delete and Find) and two more order opperations -
																																						  *             Find Next (Succ) and Find Previous (Pred). All these operations
																																						  *             are performed in O(log m) time, or equivalently in O(log log M)
																																						  *             time, where M = 2^m is the maximum number of elements that can be
																																						  *             stored in the tree. It also supports the operation Minimum (Min)
																																						  *             and Maximum (Max), wich returns the minimum and maximum element
																																						  *             stored in the tree respectively. These both run in O(1) time,
																																						  *             since the minimum and maximum element are stored as attributes in
																																						  *             each tree.
																																						  ******************************************************************************/
struct TvEB
{
	/*************************************************************************//**
																			   * @brief      Constructor.
																			   *
																			   * @param[in]  uniSize  The size of the tree universe
																			   ****************************************************************************/
	TvEB(int64_t uniSize);

	/*************************************************************************//**
																			   * @brief      Destructor.
																			   ****************************************************************************/
	~TvEB();

	/*************************************************************************//**
																			   * @brief      The size of the universe.
																			   ****************************************************************************/
	const int64_t uni;

	/*************************************************************************//**
																			   * @brief      The square root of the universe size.
																			   ****************************************************************************/
	const int64_t uniSqrt;

	/*************************************************************************//**
																			   * @brief      The lower square root of the universe size.
																			   ****************************************************************************/
	const int64_t lowerUniSqrt;

	/*************************************************************************//**
																			   * @brief      The higher square root of the universe size.
																			   ****************************************************************************/
	const int64_t higherUniSqrt;

	/*************************************************************************//**
																			   * @brief      The minimal value in the tree.
																			   ****************************************************************************/
	int64_t min;

	/*************************************************************************//**
																			   * @brief      The maximal value in the tree.
																			   ****************************************************************************/
	int64_t max;

	/*************************************************************************//**
																			   * @brief      The pointer to the summary structure of the tree.
																			   ****************************************************************************/
	TvEB * summary;

	/*************************************************************************//**
																			   * @brief      The pointer to the array of clusters of the tree.
																			   ****************************************************************************/
	TvEB ** cluster;
};

/***************************************************************************//**
																			 * @brief      Rounds up the given value to the next higher power of two.
																			 *
																			 * @param[in]  val   The value to round up.
																			 *
																			 * @return     The higher power of two from the given value.
																			 ******************************************************************************/
int64_t powTwoRoundUp(int64_t val);

/***************************************************************************//**
																			 * @brief      Returns the lower square root of the given value, which is 2^(
																			 *             floor ( log_2 ( value ) / 2 ) ).
																			 *
																			 * @param[in]  val   The value from which the lower square root is calculated.
																			 *
																			 * @return     The lower square root of the value.
																			 ******************************************************************************/
double lowerSqrt(int64_t val);

/***************************************************************************//**
																			 * @brief      Returns the higher square root of the given value, which is 2^(
																			 *             ceil ( log_2 ( value ) / 2 ) ).
																			 *
																			 * @param[in]  val   The value from which the higher square root is calculated.
																			 *
																			 * @return     The higher square root of the value.
																			 ******************************************************************************/
double higherSqrt(int64_t val);

/***************************************************************************//**
																			 * @brief      Returns the element's index in the cluster.
																			 *
																			 * @param[in]  tree  The pointer to the van Emde Boas tree.
																			 * @param[in]  val   The value of the element.
																			 *
																			 * @return     The element's index in the cluster.
																			 ******************************************************************************/
int64_t low(TvEB * tree, int64_t val);

/***************************************************************************//**
																			 * @brief      Returns the index of the element's cluster.
																			 *
																			 * @param[in]  tree  The pointer to the van Emde Boas tree.
																			 * @param[in]  val   The value of the element.
																			 *
																			 * @return     The index of the element's cluster.
																			 ******************************************************************************/
int64_t high(TvEB * tree, int64_t val);

/***************************************************************************//**
																			 * @brief      Returns the value on index low in the cluster high.
																			 *
																			 * @param[in]  tree  The pointer to the van Emde Boas tree.
																			 * @param[in]  high  The cluster index.
																			 * @param[in]  low   The index of the element in the cluster.
																			 *
																			 * @return     The value of the element.
																			 ******************************************************************************/
int64_t index(TvEB * tree, int64_t high, int64_t low);

/***************************************************************************//**
																			 * @brief      Finds the lowest value stored in the given tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[out] res    The lowest element.
																			 *
																			 * @retval     true   Successfully found the minimum.
																			 * @retval     false  Failed to found the minimum.
																			 ******************************************************************************/
bool vEB_min(TvEB * tree, int64_t & res);

/***************************************************************************//**
																			 * @brief      Finds the highest value stored in the given tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[out] res    The highest element.
																			 *
																			 * @retval     true   Successfully found the maximum.
																			 * @retval     false  Failed to found the maximum.
																			 ******************************************************************************/
bool vEB_max(TvEB * tree, int64_t & res);

/***************************************************************************//**
																			 * @brief      Inserts the given value into the given vEB tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[in]  val    The value of the element to insert.
																			 *
																			 * @retval     true   Successfully inserted the value.
																			 * @retval     false  Failed to insert the value.
																			 ******************************************************************************/
bool vEB_insert(TvEB *& tree, int64_t val, int64_t parentUniSqrt = 65536);

/***************************************************************************//**
																			 * @brief      Removes the given value from the given vEB tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[in]  val    The value of the element to remove.
																			 *
																			 * @retval     true   Successfully removed the value.
																			 * @retval     false  Failed to remove the value.
																			 ******************************************************************************/
bool vEB_delete(TvEB *& tree, int64_t val);

/***************************************************************************//**
																			 * @brief      Finds if the given value is in the given vEB tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[in]  val    The value of the element to find.
																			 *
																			 * @retval     true   Successfully found the element.
																			 * @retval     false  Failed to found the element.
																			 ******************************************************************************/
bool vEB_find(TvEB * tree, int64_t val);

/***************************************************************************//**
																			 * @brief      Finds the smallest value at least the given value in the given
																			 *             tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[in]  val    The lower bound for the value of the sought element
																			 * @param[out] res    The found element.
																			 *
																			 * @retval     true   Successfully found the successor.
																			 * @retval     false  Failed to found the successor.
																			 ******************************************************************************/
bool vEB_succ(TvEB * tree, int64_t val, int64_t & res);

/***************************************************************************//**
																			 * @brief      Finds the largest value at most the given value in the given
																			 *             tree.
																			 *
																			 * @param[in]  tree   The pointer to the van Emde Boas tree.
																			 * @param[in]  val    The upper bound for the value of the sought element
																			 * @param[out] res    The found element.
																			 *
																			 * @retval     true   Successfully found the predecessor.
																			 * @retval     false  Failed to found the predecessor.
																			 ******************************************************************************/
bool vEB_pred(TvEB * tree, int64_t val, int64_t & res);

/***************************************************************************//**
																			 * @brief      Prints pointer values of the given tree.
																			 *
																			 * @param[in]  tree  The pointer to the van Emde Boas tree.
																			 * @param      os    The output stream in which the values are printed.
																			 ******************************************************************************/
void vEB_print(TvEB * tree, std::ostream & os);

#endif /* __VEB_H_458976543568798867538649758687654752463747856374562543646__ */