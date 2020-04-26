// Сбалансированные деревья поиска.
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _TREAP_Tree_BALANCED_TREE_CPP_
#define _TREAP_Tree_BALANCED_TREE_CPP_ 1

// Рандомизированное дерево двоичного поиска.
// Достоинства: Простота и понятность реализации при обеспечении логарифмической сложности операций.
// 24 august 2017
// По материалам статьи с хабрахабра: 7 июня 2012 в 23:43 Рандомизированные деревья поиска

// Дерамида (объединение дерева и двоичной кучи). Сидель (Siedel) и Арагон (Aragon) 1996г.
// Имеет опыт использования в качестве базовой структуры в тщательно проработанной библиотеке LEDA.
// Описание заимствовано с сайта GeeksForGeeks.
// Amat Treap Node
struct TreapNode
{
	integer  priority = -1;
	data_BalTree key = data_BalTree(-1, -1);
	TreapNode* left = nullptr, * right = nullptr;
};

/* T1, T2 and T3 are subtrees of the tree rooted with y
(on left side) or x (on right side)
	y                               x
   / \     Right Rotation          /  \
  x   T3   – – – – – – – >        T1   y
 / \       < - - - - - - -            / \
T1  T2     Left Rotation            T2  T3 */

// Amat utility function to right rotate subtree rooted with y
// See the diagram given above.
TreapNode* rightRotate(TreapNode* y)
{
	TreapNode* x = y->left, * T2 = x->right;

	// Perform rotation
	x->right = y;
	y->left = T2;

	// Return new root
	return x;
}

// Amat utility function to left rotate subtree rooted with x
// See the diagram given above.
TreapNode* leftRotate(TreapNode* x)
{
	TreapNode* y = x->right, * T2 = y->left;

	// Perform rotation
	y->left = x;
	x->right = T2;

	// Return new root
	return y;
}

/* Utility function to add a new key */
TreapNode* newNode(data_BalTree key)
{
	TreapNode* temp = new TreapNode;
	temp->key = key;
	//temp->priority = rand() % 100;
	//temp->priority = rand() % 200;
	//temp->priority = rand() % 65536;
	temp->priority = rand();
	temp->left = temp->right = nullptr;
	return temp;
}

// C function to search a given key in a given BST
TreapNode* search_recursive(TreapNode* root, data_BalTree key)
{
	// Base Cases: root is null or key is present at root
	if (root == nullptr || (root->key == key)) {
		return root;
	}
	else if (root->key < key) {

		// Key is greater than root's key
		return search_recursive(root->right, key);
	}
	else {
		// Key is smaller than root's key

		return search_recursive(root->left, key);
	}

}

// Итеративный вариант должен быть быстрее.
// C function to search a given key in a given BST
TreapNode* search(TreapNode* root, data_BalTree key)
{

	TreapNode* scan = root;

	for (;;) {

		// Base Cases: root is null or key is present at root
		if (scan == nullptr || (scan->key == key)) {
			return scan;
		}
		else if (scan->key < key) {

			// Key is greater than root's key

			scan = scan->right;
		}
		else {

			// Key is smaller than root's key

			scan = scan->left;
		}
	}
}


/* Recursive implementation of insertion in Treap */
TreapNode* insert(TreapNode* root, data_BalTree key)
{
	// If root is nullptr, create a new node and return it
	if (!root)
		return newNode(key);

	// If key is smaller than root
	if (key < root->key)
	{
		// Insert in left subtree
		root->left = insert(root->left, key);

		// Fix Heap property if it is violated
		if (root->left->priority > root->priority)
			root = rightRotate(root);
	}
	else // If key is greater
	{
		// Insert in right subtree
		root->right = insert(root->right, key);

		// Fix Heap property if it is violated
		if (root->right->priority > root->priority)
			root = leftRotate(root);
	}
	return root;
}

/* Recursive implementation of Delete() */
TreapNode* deleteNode(TreapNode* root, data_BalTree key)
{
	if (root == nullptr)
		return root;

	if (key < root->key)
		root->left = deleteNode(root->left, key);
	else if (key > root->key)
		root->right = deleteNode(root->right, key);

	// IF KEY IS AT ROOT

	// If left is nullptr
	else if (root->left == nullptr)
	{
		TreapNode* temp = root->right;
		delete(root);
		root = temp;  // Make right child as root
	}

	// If Right is nullptr
	else if (root->right == nullptr)
	{
		TreapNode* temp = root->left;
		delete(root);
		root = temp;  // Make left child as root
	}

	// If ksy is at root and both left and right are not nullptr
	else if (root->left->priority < root->right->priority)
	{
		root = leftRotate(root);
		root->left = deleteNode(root->left, key);
	}
	else
	{
		root = rightRotate(root);
		root->right = deleteNode(root->right, key);
	}

	return root;
}

// Amat utility function to print tree
void inorder(TreapNode* root)
{
	if (root)
	{
		inorder(root->left);
		printf("key: i= %lld, count_neighbour = %lld  | priority: %lld ", root->key.i, root->key.count_neighbour, root->priority);
		if (root->left)
			printf(" | left child: i= %lld, count_neighbour = %lld ", root->left->key.i, root->key.count_neighbour);
		if (root->right)
			printf(" | right child: i= %lld, count_neighbour = %lld \n", root->right->key.i, root->key.count_neighbour);
		inorder(root->right);
	}
}


// Полное удаление бинарного радномизированного дерева поиска.
void clear_random_tree(TreapNode* root)
{
	if (root != nullptr) {
		clear_random_tree(root->left);
		clear_random_tree(root->right);
		// удаляем лист.
		delete root;
		//free(p);
		root = nullptr;
	}
} // clear_random_tree

TreapNode* findmax_random_tree(TreapNode*& p)
{
	TreapNode* p1 = p;
	// поиск узла с минимальным ключом в дереве p
	if (p1 != nullptr) {
		//return p1->right ? findmax_random_tree(p1->right): p1;
		while (p1->right != nullptr) p1 = p1->right;
		TreapNode* q = new TreapNode;
		q->key = p1->key;
		p1 = nullptr;
		return q;
	}
	else {
		// На поиск максимума подан нулевой указатель.
		return nullptr;
	}
} // findmax

// Driver Program to test above functions
/*
int test_Treap()
{
	srand(time(nullptr));

	struct TreapNode *root = nullptr;
	root = insert(root, 50);
	root = insert(root, 30);
	root = insert(root, 20);
	root = insert(root, 40);
	root = insert(root, 70);
	root = insert(root, 60);
	root = insert(root, 80);

	printf("Inorder traversal of the given tree \n");
	inorder(root);

	printf("\nDelete 20\n";
	root = deleteNode(root, 20);
	printf("Inorder traversal of the modified tree \n");
	inorder(root);

	printf("\nDelete 30\n";
	root = deleteNode(root, 30);
	printf("Inorder traversal of the modified tree \n");
	inorder(root);

	printf("\nDelete 50\n");
	root = deleteNode(root, 50);
	printf("Inorder traversal of the modified tree \n");
	inorder(root);

	TreapNode *res = search(root, 50);
	(res == nullptr) ? printf("\n50 Not Found "):
		printf("\n50 found");

	return 0;
}
*/

/*
struct node_random_tree // структура для представления узлов дерева
{//+
	data_BalTree key;
	integer size;
	node_random_tree* left;
	node_random_tree* right;
	node_random_tree(data_BalTree k) { key = k; left = right = 0; size = 1; }
};

//+
// Обычный двоичный поиск.
node_random_tree* find_random_tree(node_random_tree* p, data_BalTree k) // поиск ключа k в дереве p
{


	node_random_tree* q = p;

	for (;;) {
		if (p == nullptr) break;  // не найден.
		else if (k < p->key) {
			p = p->left;
		}
		else if (k > p->key) {
			p = p->right;
		}
		else break; // найден.
	}

	if (p == nullptr) {
		p = q;
		q = nullptr;
		return nullptr;
	}
	else {
		node_random_tree* q1 = new node_random_tree(p->key);
		p = q;
		q = nullptr;
		return q1;
	}

} // find_random_tree

//+
integer getsize_random_tree(node_random_tree* p) // обертка для поля size, работает с пустыми деревьями (t=nullptr)
{
	if (p==nullptr) return 0;
	return p->size;
}

//+
void fixsize_random_tree(node_random_tree* &p) // установление корректного размера дерева
{
	if (p!=nullptr) p->size = getsize_random_tree(p->left) + getsize_random_tree(p->right) + 1;
}

//+
void insert_random_tree(node_random_tree* &p, data_BalTree k) // классическая вставка нового узла с ключом k в дерево p
{
	if (p == nullptr) {
		p = new node_random_tree(k);
	}
	else {

		if (k < p->key)
			 insert_random_tree(p->left, k);
		else
			 insert_random_tree(p->right, k);

		fixsize_random_tree(p);

	}
} // insert_random_tree

//+
void rotateright_random_tree(node_random_tree* &p) // правый поворот вокруг узла p
{
	if (p == nullptr || p->left == nullptr) {
		// ничего не делаем.
	}
	else {
		node_random_tree* q = p->left;
		if (q != nullptr) {
			p->left = q->right;
			q->right = p;
			fixsize_random_tree(p);
			fixsize_random_tree(q);
		}
	}
} // rotateright_random_tree

//+
void rotateleft_random_tree(node_random_tree* &q) // левый поворот вокруг узла q
{
	if (q == nullptr || q->right == nullptr) {
		// ничего не делаем.
	}
	else {
		node_random_tree* p = q->right;
		if (p != nullptr) {
			q->right = p->left;
			p->left = q;
			p->size = q->size;
			fixsize_random_tree(p);
			fixsize_random_tree(q);
		}
	}
} // rotateleft_random_tree

void insertroot_random_tree(node_random_tree* &p, data_BalTree k) // вставка нового узла с ключом k в корень дерева p
{
	if (p == nullptr) {
		//return new node_random_tree(k);
		p=new node_random_tree(k);
	}
	else {

		if (k < p->key) {
			insertroot_random_tree(p->left, k);
			rotateright_random_tree(p);
		}
		else  {
			insertroot_random_tree(p->right, k);
			rotateleft_random_tree(p);
		}

		fixsize_random_tree(p);
	}

}// insertroot_random_tree


void insert_work_random_tree(node_random_tree* &p, data_BalTree k) // рандомизированная вставка нового узла с ключом k в дерево p
{
	if (p == nullptr) { p = new node_random_tree(k);  }
	else {
		if (rand() < RAND_MAX / (p->size + 1)) {
			insertroot_random_tree(p, k);
		}
		else {

			if (k < p->key)
				insert_work_random_tree(p->left, k);
			else if (k > p->key)
				insert_work_random_tree(p->right, k);

		}
		fixsize_random_tree(p);
	}
} // insert_work_random_tree

node_random_tree* join_random_tree_base(node_random_tree* &p, node_random_tree* &q) // объединение двух деревьев
{
	if (p == nullptr) return q;
	if (q == nullptr) return p;

	insert_work_random_tree(q, p->key);
	q->left = join_random_tree_base(p->left,q->left);
	q->right= join_random_tree_base(p->right, q->right);
	p->left = nullptr;
	p->right = nullptr;
	delete p;
	p = nullptr;
	fixsize_random_tree(q);
	return q;
}

node_random_tree* join_random_tree(node_random_tree* &p, node_random_tree* &q) // объединение двух деревьев
{
	if (p==nullptr) return q;
	if (q==nullptr) return p;
	if (rand() /(RAND_MAX/(p->size + q->size) +1)<p->size)
	{
		p->right = join_random_tree_base(p->right, q);
		fixsize_random_tree(p);
		return p;
	}
	else
	{
		q->left = join_random_tree_base(p, q->left);
		fixsize_random_tree(q);
		return q;
	}
} // join_random_tree
*/
/*
node_random_tree* remove_random_tree(node_random_tree* p, data_BalTree k) // удаление из дерева p первого найденного узла с ключом k
{
	if (p==nullptr) return p;

	if (k < p->key) {
		p->left = remove_random_tree(p->left, k);
	}
	else if (k > p->key) {
		p->right = remove_random_tree(p->right, k);
	}
	else {
		// ключи равны.
		node_random_tree* q = join_random_tree(p->left, p->right);
		delete p;
		return q;
	}

	return p;
}
*/
/*
void remove_random_tree(node_random_tree* &p, data_BalTree k) //удаление из дерева p первого найденного узла с ключом k
{
	if (p != nullptr) {


		node_random_tree* p1 = p;
		//node_random_tree* p1_parent = p;
		bool bleft = true;


		for (;;) {
			if (p1 != nullptr) {
				if (k < p1->key) {
					bleft = true;
					//p1_parent = p1;
					 p1 = p1->left;
				}
				else if (k > p1->key) {
					bleft = false;
					//p1_parent = p1;
					p1 = p1->right;
				}
				else if (p1->key == k)  {

				/*if (p1_parent == p1) {

					if (p1->key == k) {
						p1 = nullptr;
						p1_parent = nullptr;
						delete p;
						p = nullptr;
					}
				}
				else {*/
				/*
				// найден.
				// ключи равны.
				node_random_tree* p1_mem = p1;
				//node_random_tree* q = join_random_tree(p1->left, p1->right);
				p1= join_random_tree(p1->left, p1->right);
				*/
				/*
				if (bleft) {
					p1_parent->left = q;
					q = nullptr;
				}
				else {
					p1_parent->right = q;
					q = nullptr;
				}
				p1_mem->left = nullptr;
				p1_mem->right = nullptr;
				delete p1_mem;
				*/
				/*
									p1_mem->left = nullptr;
									p1_mem->right = nullptr;
									delete p1_mem;
									p1_mem = nullptr;
									//p1_parent = nullptr;
									//p1 = q;
									// досрочный выход из цикла for.
									break;
								//}
							}
							}
							else {
								// досрочный выход из цикла for.
								break;
							}
						}

					}
				}

				// Полное удаление бинарного радномизированного дерева поиска.
				void clear_random_tree(node_random_tree* &p)
				{
					if (p != nullptr) {
						clear_random_tree(p->left);
						clear_random_tree(p->right);
						// удаляем лист.
						delete p;
						//free(p);
						p = nullptr;
					}
				} // clear_random_tree

				node_random_tree* findmax_random_tree(node_random_tree* &p)
				{
					node_random_tree* p1 = p;
					// поиск узла с минимальным ключом в дереве p
					if (p1 != nullptr) {
						//return p1->right ? findmax_random_tree(p1->right): p1;
						while (p1->right != nullptr) p1 = p1->right;
						node_random_tree* q = new node_random_tree(p1->key);
						p1 = nullptr;
						return q;
					}
					else {
						// На поиск максимума подан нулевой указатель.
						return nullptr;
					}
				} // findmax
				*/
				// Конец данных рандомизированное дерево двоичного поиска.

#endif