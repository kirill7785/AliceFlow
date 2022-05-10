#pragma once
#ifndef _MY_HANDLE_ERROR_CPP_
#define _MY_HANDLE_ERROR_CPP_ 1

/*
	Проверяет работу функции handle_error.
В случае если malloc возвращает nullptr то 
еще дважды программа пытается выделить  данный объем памяти заданного типа.
В случае неудачи осуществляется выход из программы. 
Перед каждым новым выделением памяти и перед выходом из программы организуется пауза
в работе программы для пользователя.

Предполагается что выделяемый объем памяти хорошо известен заранее, выделяется
ровно столько сколько нужно и меньшее значение объёма выделяемой памяти никак не подойдет.

04.09.2020

*/



/*
void handle_error(Ak1* &handle, char* ch_var, char* ch_function, integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (Ak1*)malloc((n) * sizeof(Ak1));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (Ak1*)malloc((n) * sizeof(Ak1));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}

void handle_error(Ak1* &handle, const char ch_var[], const char ch_function[], integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (Ak1*)malloc((n) * sizeof(Ak1));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (Ak1*)malloc((n) * sizeof(Ak1));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}

void handle_error(Ak* &handle, char* ch_var, char* ch_function, integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (Ak*)malloc((n) * sizeof(Ak));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (Ak*)malloc((n) * sizeof(Ak));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}


void handle_error(Ak* &handle, const char ch_var[], const char ch_function[], integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (Ak*)malloc((n) * sizeof(Ak));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (Ak*)malloc((n) * sizeof(Ak));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}
*/

template <typename doublerealT>
void handle_error(doublerealT*& handle, char* ch_var, char* ch_function, integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}


template <typename doublerealT>
void handle_error(doublerealT*& handle, const char ch_var[], const char ch_function[], integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}

template <typename doublerealT>
void handle_error(doublerealT*& handle, char* ch_var, integer id_var, char* ch_var2, char* ch_function, integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
		printf("Problem: not enough memory on your equipment for %s %lld %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#else
		printf("Problem: not enough memory on your equipment for %s %d %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#endif
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
		if (handle == nullptr) {
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for %s %lld %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#else
			printf("Problem: not enough memory on your equipment for %s %d %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#endif
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
			if (handle == nullptr) {
#if doubleintprecision == 1
				printf("Problem: not enough memory on your equipment for %s %lld %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#else
				printf("Problem: not enough memory on your equipment for %s %d %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#endif
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
#if doubleintprecision == 1
				printf("ok! memory allocation successfully for %s %lld %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#else
				printf("ok! memory allocation successfully for %s %d %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#endif
			}
		}
		else {
#if doubleintprecision == 1
			printf("ok! memory allocation successfully for %s %lld %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#else
			printf("ok! memory allocation successfully for %s %d %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#endif
		}
	}
}

template <typename doublerealT>
void handle_error(doublerealT*& handle, const char ch_var[], integer id_var, const char ch_var2[], const char ch_function[], integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
#if doubleintprecision == 1
		printf("Problem: not enough memory on your equipment for %s %lld %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#else
		printf("Problem: not enough memory on your equipment for %s %d %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#endif
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
		if (handle == nullptr) {
#if doubleintprecision == 1
			printf("Problem: not enough memory on your equipment for %s %lld %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#else
			printf("Problem: not enough memory on your equipment for %s %d %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#endif
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (doublerealT*)malloc((n) * sizeof(doublerealT));
			if (handle == nullptr) {
#if doubleintprecision == 1
				printf("Problem: not enough memory on your equipment for %s %lld %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#else
				printf("Problem: not enough memory on your equipment for %s %d %s in %s ...\n", ch_var, id_var, ch_var2, ch_function);
#endif
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
#if doubleintprecision == 1
				printf("ok! memory allocation successfully for %s %lld %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#else
				printf("ok! memory allocation successfully for %s %d %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#endif
			}
		}
		else {
#if doubleintprecision == 1
			printf("ok! memory allocation successfully for %s %lld %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#else
			printf("ok! memory allocation successfully for %s %d %s in %s.\n", ch_var, id_var, ch_var2, ch_function);
#endif
		}
	}
}

/*
void handle_error(bool* handle, char* ch_var, char* ch_function, integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (bool*)malloc((n) * sizeof(bool));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (bool*)malloc((n) * sizeof(bool));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}


void handle_error(bool* handle, const char ch_var[], const char ch_function[], integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (bool*)malloc((n) * sizeof(bool));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (bool*)malloc((n) * sizeof(bool));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}

void handle_error(integer* handle, char* ch_var, char* ch_function, integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (integer*)malloc((n) * sizeof(integer));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (integer*)malloc((n) * sizeof(integer));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}

void handle_error(integer* handle, const char ch_var[], const char ch_function[], integer n)
{
	if (handle == nullptr) {
		// недостаточно памяти на данном оборудовании.
		printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
		printf("Please any key to continue...\n");
		system("pause");
		// Здесь пользователь может закрыть другие приложения и освободить память.
		// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
		handle = (integer*)malloc((n) * sizeof(integer));
		if (handle == nullptr) {
			printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
			printf("Please any key to continue...\n");
			system("pause");
			// Здесь пользователь может закрыть другие приложения и освободить память.
			// После этого память станет доступной и программа сможет в теории продолжить своё выполнение.
			handle = (integer*)malloc((n) * sizeof(integer));
			if (handle == nullptr) {
				printf("Problem: not enough memory on your equipment for %s in %s ...\n", ch_var, ch_function);
				printf("Please any key to exit...\n");
				system("pause");
				exit(1);
			}
			else {
				printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
			}
		}
		else {
			printf("ok! memory allocation successfully for %s in %s.\n", ch_var, ch_function);
		}
	}
}
*/

#endif