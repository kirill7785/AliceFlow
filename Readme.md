# AliceFlow_v0.48 - Программа для тепловых расчётов в 3D.

## Применение
Модели физических тел с внутренними и внешними источниками тепла и различной организацией теплоотводов. Возможен расчёт систем охлаждения с воздухом и жидкостью.

Программа Alice_Flow_v0.48 предназначена для расчёта поля температур в трёхмерных твердотельных моделях. В ряде случаев учитывается конвективный перенос теплоносителя. Учитываются нелинейности  различной природы, например, зависимость теплопроводности материала от температуры. Поддерживается расчёт тепловой переходной характеристики - зависимости максимальной температуры в расчётной области от времени. Для ускорения вычислений используется алгебраический многосеточный метод. Для сокращения числа ячеек сетки и ускорения нестационарных  расчётов реализованы Адаптивные Локально Измельчённые Сетки (АЛИС).

Программа AliceFlow.v.0.48 решает систему уравнений в частных производных. Для расчета пользователю необходимо задать: 1. Геометрические примитивы из которых состоит тепловая модель. Поддерживаются кубики с прямыми углами, цилиндры, полигоны а также CAD объекты из stl бинарных файлов. Примитивы располагаются в дереве модели. Если два примитива пересекаются то область пересечения принадлежит примитиву расположенному в дереве позднее (более высокий приоритет). 
2. Свойства материалов (постоянные или зависящие от температуры):
* для стационарной теплопередачи в твёрдом теле только коэффициент теплопроводности;
* для нестационарной теплопередачи в твердом теле  коэффициент теплопроводности, плотность, удельную теплоёмкость при постоянном давлении;
* Для расчёта течения изотермической вязкой жидкости требуется задать плотность и динамическую вязкость;
* Для расчёта течения неизотермической жидкости в приближении Обербека -Буссинеска требуется задать плотность, динамическую вязкость, теплопроводность жидкости, удельную теплоёмкость при постоянном давлении в жидкости, коэффициент линейного теплового расширения в жидкости. Для решения задач сопряженного теплоообмена в кажом твёрдом теле требуется задать плотность, коэффициент теплопроводности, удельную теплоёмкость при постоянном давлении. 
* Для расчёта стационарного напряженно деформируемого состояния твердого тела требуется задать модуль Юнга и коэффициент Пуассона.
3. Тепловые мощности в Вт в некоторых геометрических примитивах по желанию пользователя. Можно задать постоянную мощность не зависящую от времени, а можно зависящую от времени в соответствии с предопределенными законами. Поддерживается задание импульсных режимов (square wave) - чередование включений и выключений мощности, а также пользовательских законов изменения мощности от времени piecewise const. 
4. (граичные условия). Теплоотводы, например, стенки (wall) с заданной температурой. При гидродинамическом расчёте необходимо задать скорость и температуру на входе, границы с условиями симметрии и выходные границы.
  
 ## Ссылка на Realese (исполняемый файл)
Для пользователей которые не хотят лезть в исходные тексты программы AliceFlow и заниматься их компиляцией для ОС Windows x64 приготовлен скомпилированный исполняемый файл (три исполняемых exe файла) которым пользуется автор программы. Программа имеет интуитивно понятный графический тньнофейс пользователя в стиле пограммы ANSYS Icepak. С программой AliceFlow предоставляется набор простейших моделей на которых сразу можно испытать программу в действии.
Ссылка на realese:
https://github.com/kirill7785/AliceFlow/tree/master/Alice%20EXE
Запускать нужно exe файл графического интерфейса пользователя ProjectLaplas.exe. Сам солвер AliceFlow_v0_48.exe должен при этом находится в папке ./test_pattern/solver/x64/AliceFlow_v0_48.exe
  
## Реализованные алгоритмы
 
* 3D Температурный решатель как в твёрдом теле так и для задач сопряженного теплообмена. Метод контрольного объёма.
Стационарный или нестационарный решатели. Для задач с конвекцией поле скорости можно задать аналитически, либо можно
загрузить из текстового файла load.txt, предварительно рассчитав его на основе встроенного SIMPLE алгоритма. Также
поле скорости в модели можно найти на основе решения уравнения Эйлерова потенциального течения (inviscid flow model).
Для задач расчёта температуры в твёрдом теле поддерживаются также граничные условия Ньютна-Рихмана (обычно если поле скорости не задано)
или Стефана-Больцмана (для теплообмена излучением).
* 3D cfd Полунеявная процедура для связывания уравнений неразрывности и скорость-давление  (SIMPLE [1972]). Стационарный или неастационарный гидродинамический решатель доступны.
* Поправка С.М. Рхи и У.Л. Чоу [1984] См. статью Ибрагима Сезая: https://github.com/kirill7785/AliceFlow/blob/master/I.%20Sezai%20перевод%20с%20английского.pdf
* Адаптивные Локально Измельченные Сетки (АЛИС) (неструктурированная сетка). Из меню доступны настройки подробности сетки.
Адаптивная сетка генерируется автоматическим образом без участия пользователя (автомат).
* Схемы высокой разрешающей способности для аппроксимации конвекции как на структурированной так и на АЛИС неравномерной сетке: WACEB, SMARTER, SUPER-C и др. 
* Алгебраический многосеточный метод Джона Руге и Клауса Штубена. BiCGStab. ilu0 сглаживатель и др.
Поддерживаются две библиотеки : CUSP NVIDIA 0.5.1 или AMGCL ddemidov реализующие параллельные версии
алгебраического многосеточного метода сглаженной агрегации усиленные BiCGStab или FGMRes(m) алгоритмами.
Поддерживается Си версия алгоритма amg1r5 Джона Руге и Клауса Штубена, amg1r6 версия также, + BiCGStab 
или FGMRes(m) алгоритмы для amg1r5(r6), + ilu(k) сглаживатель на основе двочной кучи на массиве для amg1r5(r6).
Также реализован собственный алгебраический многосеточный метод РУМБАv.0.14 на основе RS, PMIS или HMIS огрубления.
Алгоритм PMIS полностью распараллелен. Поддерживается усиление многосеточного метода BiCGStab или FGMRes(m) алгоритмами, +ilu(k) сглаживатель.
Используется усечение оператора интерполяции, Порог сильной связи зависящий от номера уровня, смешанная точность
вещественной арифметики float для V -цикла и double для BiCGStab или FGMRes(m). Для RS огрубления реализована 
быстродействующая приоритетная очередь на Фибоначчиевой куче.
* Модели турбулентности : Спаларт Аллмарес, SST K-Omega Ментер.
* Модель ламинарно-турбулентного перехода Ментора и Латгрии [2009] gamma-ReTheta-SST model. (Transition SST model 4 eqn).
* Приближение Обербека-Буссинеска. Сопряженный теплообмен при естественно конвективном охлаждении.
* Криволинейная граница расчётной области аппроксимируется ступеньками. Метод подсеточного разрешения геометрии не реализован (отсутствует). 
Рекомендуется использовать прямоугольные расчётные области.
* Решатель написан на С/C++ и поддерживает распараллеливание на 8 потоков с помощью технологии OpenMP.
* Дружелюбный графический интерфейс пользователя на  Embarcadero® Delphi XE8 Version 22.0.19027.8951.
* Графовый метод решения уравнения теплопроводности в нестационарной постановке. Поддерживается нелинейное граничное условие Стефана-Больцмана и 
модель зазора (зазоров) с теплообменом излучением между границами зазора. СЛАУ решается двумя алгебраическими многосеточными методами на выбор: amg1r5 или РУМБА.v.0.14.
* Графический визуализатор на GLFW OpenGL на языке с++, встроенный в код программы. Визуализация результатов расчёта в 3D, анимация результатов нестационарного CFD расчёта. 
  
  # Системные требования для компиляции програмы из исходных текстов: 
 1. Вариант а).
* 1.1. OС Windows x64
* 1.2. компилятор: Visual Studio 2015 community
* 1.3. nvidia cuda toolkit 8.0
* 1.4. nvidia cusp library 0.5.1
* 1.5. компилировать с опцией /bigobj
* 1.6. openmp выключить. 
* 1.7. библиотека GLFW OpenGL
 2. Вариант б). 
* 2.1. OС Windows x64
* 2.2. компилятор: Visual Studio 2017(or 2019) community
* 2.3. c++ boost 1.7.0 библиотека
* 2.4. c++ amgcl 10.01.2021 библиотека
* 2.5. компилировать с опцией /bigobj
* 2.6. openmp выключить или включить 
* 2.7. библиотека GLFW OpenGL
 3. Вариант в). 
* 3.1. Программа собирается компилятором GNU g++ (g++ 9.1). C:\AliceFlow_v0_48>g++ AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt 
* 3.2. 04.08.2019 с подключенной библиотекой amgcl Дениса Демидова.
* 3.3. библиотека GLFW OpenGL
   
4. Для работоспособности exe программы консольного солвера на компьютерах под управлением ОС Windows без установленной Visual Studio необходимо скачать и установить 64 битную версию -microsoft redistributable package x64 VC_redist.x64.exe
   
5. Для визуализации результатов вычисления необходимо установить
https://www.tecplot.com/products/tecplot-360/ 
 или 
https://www.paraview.org/download/
  
## Быстрая инструкция по сборке проекта из исходных текстов

(описание устарело. Проект собирается в Microsoft Visual Studio 2019. Библиотека GLFW OpenGL должна быть прописана в зависимостях.).

### Windows

1. Вы работаете на компьютере под управлением ОС Windows 10 от Microsoft. Далее рассматривается использование только
свободных (opensource) программных средств.
2. Установите компилятор MinGW GNU g++. Для этого нужно подключение к сети интернет на вашем компьютере.
Как установить компилятор MinGW (GCC/G++) Compiller in Windows 10
https://www.youtube.com/watch?v=sXW2VLrQ3Bs
3. Скачайте с GitHub https://github.com/kirill7785/AliceFlow 
папку src и скопируйте её на диск С. Перейдите в папку C:\src.
4. Запустите Windows PowerShell. Перейдите в папку с исходным кодом cd C:\src.
5. Введите терминале g++ -o AliceFlow_v0_48.exe AliceFlow_v0_48.cpp 2> gcc_log.txt В результате появится файл AliceFlow_v0_48.exe.
6. Создайте на рабочем столе папку Alice_EXE. Поместите в папку Alice_EXE программу интерфейс AliceMesh_v0_45.exe написанную на Delphi 
(скачивается с GitHub https://github.com/kirill7785/AliceFlow).
Внутри папки Alice_EXE создайте иерархию папок test_pattern\solver\x64. В папку Alice_EXE\test_pattern\solver\x64 положите исполняемый 
файл AliceFlow_v0_48.exe. 
7. На сайте GitHub в папке AliceEXE лежат примеры для программы AliceMesh_v0_45.exe. Поместите их в папку AliceEXE на рабочем столе.
8. Можно пользоваться. Запустите AliceMesh_v0_45.exe прочитайте один из примеров File->Read. Запустите пример Solve->Run.
Программа интерфейс AliceMesh_v0_45.exe автоматически вызовет программу решатель AliceFlow_v0_48.exe.
9. Установите программу tecplot360 или бесплатный аналог - программу ParaView. Откройте в tecplot360(или ParaView) файл с расширением *.PLT
который был записан на диск после окончания работы программы AliceFlow_v0_48.exe. Вставьте картинки в отчёт в программе MS Word.
10. Для тех кто хочет использовать параллельную версию программы запустите в терминале mingw-get (Для этого нужно подключение к сети интернет
на вашем компьютере) и выберите в нём пакеты  mingw32-pthreads-w32.
11. После установки находясь в папке C:\src выполните g++ -o AliceFlow_v0_48.exe AliceFlow_v0_48.cpp -fopenmp 2> gcc_log.txt для создания 
параллельной версии программы AliceFlow_v0_48.exe.

### Linux

1. Вы работаете на компьютере под управлением ОС Linux Ubuntu-20.04.1-desktop-amdx64. http://releases.ubuntu.com/20.04/
2. Программа солвер AliceFlow_v0_48 может быть откомпилирована и работать под ОС Linux. Программа интерфес AliceMesh_v0_45.exe
работает только под ОС Windows т.к. написана на Delphi.
3. скопируйте в домашнюю папку ~/ папку src с GitHub https://github.com/kirill7785/AliceFlow.
4. Перейдите в папку src: cd ~/src.  Введите терминале g++ -o AliceFlow_v0_48 AliceFlow_v0_48.cpp 2> gcc_log.txt 
В результате в папке ~/src появится исполняемый файл AliceFlow_v0_48.
5. Если нету g++ установите его как подсказывает ОС Linux Ubuntu sudo apt get ...
6. Присвойте права на исполнение файлу AliceFlow_v0_48: chmod +x AliceFlow_v0_48.
7. Пропишите путь к исполняемому файлу в переменной PATH с помощью команды
export PATH="$HOME/src:$PATH".
8. В папке src лежит пример файла premeshin.txt - это входной файл для солвера AliceFlow_v0_48 который генерируется интерфейсом AliceMesh_v0_45.exe.
9. Вылоните в терминале находясь в папке ~/src команду AliceFlow_v0_48
10. Солвер запустится и успешно отработает на файле premeshin.txt. Для создания файлов premeshin.txt под ваши задачи нужно использовать
программу интерфейс AliceMesh_v0_45.exe работающую под управлением ОС Windows. Можно использовать, например, VMWare. Программа интерфейс 
AliceMesh_v0_45.exe не требовательна к ресурсам ПК т.к. используется только для создания модели для расчёта.
  
## Рассмотрим примеры решения задач в данной программе:

Модуль водяного охлаждения

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/water_cooling_module.png)


![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/speed.jpg)

## Адаптивная Локально Измельченная Сетка (АЛИС)

Начиная с лета 2016 года в данной программе запрограммированы Адаптивные Локально Измельченные  Сетки (АЛИС). 
При построении сетки в данной программе,  в расчетной области мы считаем, что она состоит из прямоугольных параллелепипедов с возможным локальным измельчением, диктуемым условиями задачи. 
Сущность технологии построения локально-адаптивных сеток заключается в следующем. Начальная сетка является декартовой, и все ее ячейки являются прямоугольными параллелепипедами. Затем в соответствии с заданными критериями выделяются подобласти с особенностями геометрии или решения, и в этих подобластях строится более мелкая, по сравнению с исходной, сетка. Мы считаем для определенности, что выделяемая особенность задается какой- либо поверхностью. Если расчетная ячейка лежит в зоне влияния выделенной особенности, например, пересекается поверхностью, то такая ячейка делится на 8 равных ячеек. Далее, если необходимо, ячейки делятся еще раз, и так до достижения необходимой точности. Криволинейная граница аппроксимируется ступеньками. Ячейки начальной сетки называются ячейками уровня 0, ячейки, получаемые измельчением уровня 0, называются ячейками уровня 1 и т.д.
При генерации сеток необходимо накладывать дополнительное ограничение: в окрестности каждой ячейки не должно быть ячеек, которые отличаются от неё по размерам более чем в два раза. 

[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET3.png)](https://yadi.sk/i/Fd9L_d3bAiLD7w)

Обтекание куба на основе схем высокой разрешающей способности

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Cube%20Flow.png)

Обтекание крыла Re = 3164, Pr = 0.7.

![alt_text](https://github.com/kirill7785/AliceFlow/blob/master/picture/speed_around_wing.png)

Моделирование естественной конвекции в лабораторных условиях
Ra=6.4E+7; Pr=0.7; L/H=6.

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Raley_Benar%20Natural%20Convection.png)

Температурный пограничный слой на пластине (Блазиус 1908)
 Re=9736; Pr=0.7.

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Blasius%201908.png)


Вычисление теплового сопротивления полупроводниковых структур

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET1.png)

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picFET2.png)

Вычисление теплового сопротивления диода


[![Watch the video](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/Diod.png)](https://yadi.sk/i/3AicztJ3sbb97Q)

## Классический алгебраический многосеточный метод (CAMG)

Скорость сходимости многих методов обращения матриц может быть значительно повышена за счет использования методики под названием «многосеточный метод». Многосеточный метод включает в себя проведение ранних итераций на мелкой сетке и более поздних итерации на вложенных более грубых «виртуальных» сетках. Результаты затем передаются обратно из грубых сеток на первоначальные более подробные мелкие  сетки. С числовой точки зрения, многосеточный подход предлагает значительное преимущество. Для заданной сетки конкретного конечного размера, итерационные методы являются эффективными только при снижении ошибок, которые имеют длину волны порядка шага сетки (ребро тетраэдра при тетра мешировании, ребро гексаэдра при хекса доминант мешировании). Таким образом, в то время как более короткие длины волны ошибки исчезают довольно быстро, ошибки с большей длиной волны, порядка размера расчётной области, исчезают очень медленно (катострофически медленно). Метод Многосеточный обходит эту проблему, используя ряд грубых сеток (они являются вложенными в первоначальную подробную сетку как матрёшки) таким образом, что имеющиеся компоненты вектора ошибки с большой длиной волны являются коротковолновыми (легко падавляемыми обычным итерационным методом - сглаживателем) на грубых сетках из иерархии вложенности. Для того, чтобы избежать необходимости в грубосеточном мешировании геометрии с использованием ряда различных шагов сетки, данная программа использует Алгебраический многосеточный метод.

https://github.com/kirill7785/algebraic-multigrid-v.0.14/blob/main/README.ru.md

![alt text](https://github.com/kirill7785/algebraic-multigrid-method/blob/master/picture/picPaper.png)

# собственный визуализатор на fwgl OpenGL (render)

![alt text](https://github.com/kirill7785/AliceFlow/blob/master/picture/IMG-20201218-WA0000.jpeg)
# математическое решение задач поиска напряжённо деформированного состояния в программе AliceFlow v.0.48
![alt text](https://github.com/kirill7785/AliceFlow/blob/master/picture/Деформация%20в%20задаче%20Фламана.png)
Деформации в задаче Фламана
![alt text](https://github.com/kirill7785/AliceFlow/blob/master/picture/Напряжение%20по%20Мизесу%20в%20задаче%20Фламана.png)
Напряжения по фон Мизесу в задаче Фламана
![alt text](https://github.com/kirill7785/AliceFlow/blob/master/picture/Стальная%20балка%20малая%20деформация.png)
деформация стальной балки
![alt text](https://github.com/kirill7785/AliceFlow/blob/master/picture/Стальная%20балка%20напряжение%20по%20Мизесу.png)
напряжение по Мизесу в стальной балке протяженностью 0.5м

