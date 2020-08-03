// Сбалансированные деревья поиска.
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _SPLAY_Tree_BALANCED_TREE_CPP_
#define _SPLAY_Tree_BALANCED_TREE_CPP_ 1

// splay tree 23 01 2016.
// start.splay tree

/*
//
An implementation of top-down splaying
D. Sleator <sleator@cs.cmu.edu>
March 1992

"Splay trees", or "self-adjusting search trees" are a simple and
efficient data structure for storing an ordered set.  The data
structure consists of a binary tree, without parent pointers, and no
additional fields.  It allows searching, insertion, deletion,
deletemin, deletemax, splitting, joining, and many other operations,
all with amortized logarithmic performance.  Since the trees adapt to
the sequence of requests, their performance on real access patterns is
typically even better.  Splay trees are described in a number of texts
and papers [1,2,3,4,5].

The code here is adapted from simple top-down splay, at the bottom of
page 669 of [3].  It can be obtained via anonymous ftp from
spade.pc.cs.cmu.edu in directory /usr/sleator/public.

The chief modification here is that the splay operation works even if the
item being splayed is not in the tree, and even if the tree root of the
tree is nullptr.  So the line:

t = splay(i, t);

causes it to search for item with key i in the tree rooted at t.  If it's
there, it is splayed to the root.  If it isn't there, then the node put
at the root is the last one before nullptr that would have been reached in a
normal binary search for i.  (It's a neighbor of i in the tree.)  This
allows many other operations to be easily implemented, as shown below.

[1] "Fundamentals of data structures in C", Horowitz, Sahni,
and Anderson-Freed, Computer Science Press, pp 542-547.
[2] "Data Structures and Their Algorithms", Lewis and Denenberg,
Harper Collins, 1991, pp 243-251.
[3] "Self-adjusting Binary Search Trees" Sleator and Tarjan,
JACM Volume 32, No 3, July 1985, pp 652-686.
[4] "Data Structure and Algorithm Analysis", Mark Weiss,
Benjamin Cummins, 1992, pp 119-130.
[5] "Data Structures, Algorithms, and Performance", Derick Wood,
Addison-Wesley, 1993, pp 367-375.

The following code was written by Daniel Sleator, and is released
in the public domain.
*/

integer size_splay_Tree;  /* number of nodes in the Tree_splay */
/* Not actually needed for any of the operations */

typedef struct Tree_splay_node Tree_splay;
struct Tree_splay_node {
	Tree_splay* left, * right;
	data_BalTree item;

	Tree_splay_node() {
		left = nullptr;
		right = nullptr;
	}
};

Tree_splay* findmax(Tree_splay* p)
{
	// поиск узла с минимальным ключом в дереве p
	if (p != 0) {
		return p->right ? findmax(p->right) : p;
	}
	else {
		// На поиск максимума подан нулевой указатель.
		return 0;
	}
} // findmax

// Возвращает true если узел найден в SPLAY дереве
bool isfound_recursive(Tree_splay* p, data_BalTree k)
{
	if (p == 0) return false; // ненайден.
	if (k < p->item)
		return isfound_recursive(p->left, k);
	else if (k > p->item)
		return isfound_recursive(p->right, k);
	else return true; // найден.
} // isfound_recursive in splay tree.

  


Tree_splay* splay(data_BalTree i, Tree_splay* t) {
	/* Simple top down splay, not requiring i to be in the Tree_splay t.  */
	/* What it does is described above.							 */
	Tree_splay Nbuf, * l = nullptr, * r = nullptr, * y = nullptr;
	if (t == nullptr) return t;
	Nbuf.left = Nbuf.right = nullptr;
	l = r = &Nbuf;

	for (;;) {
		if (i < t->item) {
			if (t->left == nullptr) break;
			if (i < t->left->item)
			{
				y = t->left;						   /* rotate right */
				t->left = y->right;
				y->right = t;
				t = y;
				if (t->left == nullptr) break;
			}
			r->left = t;							   /* link right */
			r = t;
			t = t->left;
		}
		else if (i > t->item) {
			if (t->right == nullptr) break;
			if (i > t->right->item) {
				y = t->right;						  /* rotate left */
				t->right = y->left;
				y->left = t;
				t = y;
				if (t->right == nullptr) break;
			}
			l->right = t;							  /* link left */
			l = t;
			t = t->right;
		}
		else {
			break;
		}
	}
	l->right = t->left;								/* assemble */
	r->left = t->right;
	t->left = Nbuf.right;
	t->right = Nbuf.left;
	return t;
}

// Возвращает true если узел найден в SPLAY дереве
bool isfound(Tree_splay* p1, data_BalTree k)
{
	//if (p == 0) return false; // не найден.

	Tree_splay* p;
	p = p1;

	for (;;) {
		if (p == 0) return false; // не найден.
		else if (k < p->item) {
			p = p->left;
		}
		else if (k > p->item) {
			p = p->right;
		}
		else {
			//p1 = splay(k, p); // так медленнее
			p = 0;
			
			return true; // найден.
		}
	}
} // isfound in splay tree.

// Данный метод не используется.
/* Here is how sedgewick would have written this.		*/
/* It does the same thing.							    */
Tree_splay* sedgewickized_splay(data_BalTree i, Tree_splay* t) {
	Tree_splay Nbuf, * l = nullptr, * r = nullptr, * y = nullptr;
	if (t == nullptr) return t;
	Nbuf.left = Nbuf.right = nullptr;
	l = r = &Nbuf;

	for (;;) {
		if (i < t->item) {
			if (t->left != nullptr && i < t->left->item) {
				y = t->left; t->left = y->right; y->right = t; t = y;
			}
			if (t->left == nullptr) break;
			r->left = t; r = t; t = t->left;
		}
		else if (i > t->item) {
			if (t->right != nullptr && i > t->right->item) {
				y = t->right; t->right = y->left; y->left = t; t = y;
			}
			if (t->right == nullptr) break;
			l->right = t; l = t; t = t->right;
		}
		else break;
	}
	l->right = t->left; r->left = t->right; t->left = Nbuf.right; t->right = Nbuf.left;
	return t;
}

Tree_splay* insert(data_BalTree i, Tree_splay* t) {
	/* Insert i into the Tree_splay t, unless it's already there.	*/
	/* Return a pointer to the resulting Tree_splay.				 */
	//Tree_splay * new_splay_node=nullptr;

	Tree_splay* new_splay_node = new Tree_splay;
	//new_splay_node = (Tree_splay *)malloc(sizeof(Tree_splay));
	//if (new_splay_node == nullptr) {
		//printf("Ran out of space in SPLAY tree\n");
		//system("pause");
		//exit(1);
	//}
	new_splay_node->item = i;
	if (t == nullptr) {
		new_splay_node->left = new_splay_node->right = nullptr;
		size_splay_Tree = 1;
		return new_splay_node;
	}
	t = splay(i, t);
	if (i < t->item) {
		new_splay_node->left = t->left;
		new_splay_node->right = t;
		t->left = nullptr;
		size_splay_Tree++;
		return new_splay_node;
	}
	else if (i > t->item) {
		new_splay_node->right = t->right;
		new_splay_node->left = t;
		t->right = nullptr;
		size_splay_Tree++;
		return new_splay_node;
	}
	else { /* We get here if it's already in the Tree_splay */
		/* Don't add it again					  */
		//free(new_splay_node);
		delete new_splay_node;
		return t;
	}
}

Tree_splay* delete_splay_Tree(data_BalTree i, Tree_splay* t) {
	/* Deletes i from the Tree_splay if it's there.			   */
	/* Return a pointer to the resulting Tree_splay.			  */
	Tree_splay* x;
	if (t == nullptr) return nullptr;
	t = splay(i, t);
	if (i == t->item) {			   /* found it */
		if (t->left == nullptr) {
			x = t->right;
		}
		else {
			x = splay(i, t->left);
			x->right = t->right;
		}
		size_splay_Tree--;
		//free(t);
		delete t;
		t = nullptr;
		return x;
	}
	return t;						 /* It wasn't there */
} // delete_splay_Tree 

// Вставка ключа К в дерево если ключа k_search
// еще нет в дереве или модификация ключа К_search на К.
// Возвращает новое значение корня АВЛ дерева.
Tree_splay* insert_and_modify(Tree_splay*& p, data_BalTree k, data_BalTree k_search)
{
	if (isfound(p, k_search) == false) {
		//   узла в дереве нет.
		p = insert(k, p);
		return p;
	}
	else {
		// удаление k_search
		p = delete_splay_Tree(k_search, p); // необходимое действие
		p = insert(k, p);							 // вставка к.
		return p;
	}
} // insert_and_modify splay tree

// Полное удаление бинарного splay дерева.
void clear_SPLAY(Tree_splay* p)
{
	if (p != 0) {
		clear_SPLAY(p->left);
		clear_SPLAY(p->right);
		// удаляем лист.
		//free(p);
		delete p;
		p = 0;
	}
} // clear_SPLAY

/*
// Для целочисленного item.
void test_splay_tree() {
	// Amat sample use of these functions.  Start with the empty Tree_splay,
	// insert some stuff into it, and then delete it
	Tree_splay * root;
	integer i;
	root = nullptr;			 // the empty Tree_splay
	size_splay_Tree_splay = 0;
	for (i = 0; i < 1024; i++) {
		root = insert((541 * i) & (1023), root);
	}
	for (i = 0; i < 1024; i++) {
		root = delete_splay_Tree((541 * i) & (1023), root);
	}
	#if doubleintprecision == 1
		printf("size_splay_Tree_splay = %d\n", size_splay_Tree_splay);
	#else
		printf("size_splay_Tree_splay = %d\n", size_splay_Tree_splay);
	#endif

}
*/

// end splay tree

#endif