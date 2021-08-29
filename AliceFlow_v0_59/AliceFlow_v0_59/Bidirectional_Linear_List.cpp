// Двунаправленный линейный список.
// Используются в алгебраическом многосеточном методе.
#pragma once
#ifndef _BIDIRECTIONAL_LINEAR_LIST_CPP_
#define _BIDIRECTIONAL_LINEAR_LIST_CPP_ 1

// Двунаправленный линейный список целочисленных элементов. 
// Может сканироваться как вперёд так и назад.
typedef struct Thashlist_i {
	integer item;
	Thashlist_i* next;
	Thashlist_i* prev;

	Thashlist_i() {
		item = -1;
		next = nullptr;
		prev = nullptr;
	}
} hashlist_i;

// вставка целочисленного ключа в список
void insertion_list_i(hashlist_i*& p, integer d) {
	if (p == nullptr) {
		p = new hashlist_i;
		if (p == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for p in insertion hashlist_i my_agregat_amg...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		p->prev = nullptr;
		p->next = nullptr;
		p->item = d;
	}
	else {
		hashlist_i* scanner = p;
		hashlist_i* r = nullptr;
		r = new hashlist_i;
		if (r == nullptr) {
			// недостаточно памяти на данном оборудовании.
			printf("Problem: not enough memory on your equipment for r in insertion hashlist_i my_agregat_amg...\n");
			printf("Please any key to exit...\n");
			//getchar();
			system("pause");
			exit(1);
		}
		r->prev = nullptr;
		r->next = scanner;
		r->item = d;
		scanner->prev = r;
		p = r;
		r = nullptr;
		scanner = nullptr;
	}
} // insertion_list_i

  // очистка списка
  // Вероятно содержит ошибку потому что память не освобождается.
  // TODO 17 dec 2015
  // Проверено 1 января 2017.
void clear_hash_list_i(hashlist_i*& p) {
	if (p != nullptr) {
		hashlist_i* head = p;
		p = nullptr;

		hashlist_i* r = head;
		while (r != nullptr) {
			head = head->next;
			if (head != nullptr)  head->prev = nullptr;
			r->next = nullptr;
			r->prev = nullptr;
			delete r;
			r = head;
		}
	}
} // clear_hash_list_i

#endif