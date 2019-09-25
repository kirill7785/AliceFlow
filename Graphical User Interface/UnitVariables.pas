unit UnitVariables;
// Переменные для удобства пользователя и будущего Optimetric`a.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, Grids, StdCtrls;

type

  // типы для вычисления параметризованного выражения
  setChar=set of char; // тип множество символов
  pTop=^Top;
  Top=record
  private
     //myoperator : string[5];
     Fmyoperator : string; // знак операции
     Fvalue : Real; // значение константы
     Fleft, Fright : pTop; // указатели на левое и правое поддерево
     procedure Setmyoperator(const valuestr : string);
     procedure Setvalue(const rval : Real);
     procedure Setleft(const lval : pTop);
     procedure Setright(const rval : pTop);
  public
      property myoperator: String read Fmyoperator write Setmyoperator;
      property value : Real read Fvalue write  Setvalue;
      property left : pTop read Fleft write  Setleft;
      property right : pTop read Fright write  Setright;
  end;



  TFormVariables = class(TForm)
    GBVariables: TGroupBox;
    StringGridVariables: TStringGrid;
    BApply: TButton;
    ButtonRename: TButton;
    procedure FormCreate(Sender: TObject);
    procedure BApplyClick(Sender: TObject);
    procedure ButtonRenameClick(Sender: TObject);
  private
    { Private declarations }

      // поиск подстроки в строке
      //function Find(const S, P : String) : Integer;
      // подстановка значений переменных в строку
      // patterns - исходная неизменяемая строка.
      // workstring - преобразованная исходная строка.
      procedure my_substitutional_of_value(patterns : String; var workstring : String);
      // рекурсивная функция конструирования поддерева
      // выражения с корнем r из строки st
      procedure Constr_Tree(var r : pTop; var st : string);
      // рекурсивное вычисление значения функции
      // если key=false, то значение не существует
      function my_count(r : pTop;  var key : Boolean) : Real;
      // освобождение оперативной памяти занимаемой двоичным деревом
      procedure my_delete(var r : pTop);



  public
    { Public declarations }

      // преобразует параметризованую строку содержащую переменные в вещественное число.
      // если bOk = true то операция преобразования прошла успешно, если bOk=false
      // то операция преобразования содержит ошибку и операцию преобразования надо отменить
      function my_real_convert(s : String; var bOk : Boolean) : Real;
      procedure all_obj_project_variable_rename(sold : String; snew : String);
      // обновление размеров после изменения значения переменной
      // которая задана пользователем.
      procedure my_update_size();
  end;

var
  FormVariables: TFormVariables;

implementation
 uses
     VisualUnit, AddSourceUnit, AddVariableUnit, UnitRenameVariable;

{$R *.dfm}

procedure Top.Setmyoperator(const valuestr: string);
  begin
     Fmyoperator:=valuestr;
  end;

  procedure Top.Setvalue(const rval: Real);
  begin
      Fvalue:=rval;
  end;

  procedure Top.Setleft(const lval: pTop);
  begin
      Fleft:=lval;
  end;

  procedure Top.Setright(const rval: pTop);
  begin
      Fright:=rval;
  end;

// инициализация формы сразу после создания.
procedure TFormVariables.FormCreate(Sender: TObject);
var
    i : Integer;
begin
   // Инициализация формы редактирования для оптиметрика :
   // Cells[столбец, строка]
   FormVariables.StringGridVariables.Cells[1,0]:='$var';
   FormVariables.StringGridVariables.Cells[2,0]:='value';
   for i:=1 to FormVariables.StringGridVariables.RowCount-1 do
   begin
      FormVariables.StringGridVariables.Cells[0,i]:=IntToStr(i);
   end;
end;

procedure TFormVariables.BApplyClick(Sender: TObject);
var
     i : Integer;
     b : Boolean;
     s : String;
begin
    // Вводит переменные в проект.
    // переменные должны идти последовательно сверху вниз.
    // переменные не удовлетворяющие этому критерию исключаются из рассмотрения.
    // имя переменной обязательно должно начинаться со знака $.
    b:=true;
    i:=0;
    while (b and (i<StringGridVariables.RowCount-2)) do
    begin
       s:=StringGridVariables.Cells[1,i+1];
       if ((length(s)>=2) and (s[1]='$')) then b:=true
       else b:=false;
       if (b) then i:=i+1; // переходим к следующей переменной
    end;
    Laplas.ivar:=i; // количество переменных.
    SetLength(Laplas.parametric,Laplas.ivar);
    for i:=0 to Laplas.ivar-1 do
    begin
       Laplas.parametric[i].svar:=Trim(StringGridVariables.Cells[1,i+1]);
       Laplas.parametric[i].sval:=Trim(StringGridVariables.Cells[2,i+1]);
    end;
    // Все переменные хранятся в таблице  StringGridVariables

    // Все последующине переменные игнорируются и их надо стереть
    for i:=Laplas.ivar+1 to StringGridVariables.RowCount-1 do
    begin
       StringGridVariables.Cells[1,i]:='';
       StringGridVariables.Cells[2,i]:='';
    end;

    // обновление размеров объектов
    my_update_size();
    // прорисовка геометрии
    Laplas.ReadyPaint;
    Close();
end;

// поиск подстроки в строке
// аналог встроенной функции pos.
(*
function TFormVariables.Find(const S, P : String) : Integer;
var
  i, j : Integer;
begin
  Result := 0;
  if Length(P) > Length(S)
   then Exit;
  for i := 1 to Length(S) - Length(P) + 1 do
  for j := 1 to Length(P) do
  if P[j] <> S[i+j-1]
  then Break
  else
    if j = Length(P)
    then
     begin
       Result := i;
       Exit;
     end;
end;
*)

// подстановка значений переменных в строку
// patterns - исходная неизменяемая строка.
// workstring - преобразованная исходная строка.
procedure TFormVariables.my_substitutional_of_value(patterns : String; var workstring : String);
var
   s : String;
   i,ir,ivar,ipos, iposfirst : Integer;
   bcont, bfirst : Boolean;
begin
   s:=patterns;
   if (Laplas.ivar>0) then
   begin
      bcont:=true;
      while (bcont) do
      begin
         // инициализация
         bcont:=false;
         ipos:=0;
         bfirst:=true;
         iposfirst:=0;
         ivar:=-1; // несуществующая переменная
         for i:=0 to Laplas.ivar-1 do
         begin
            ir:=Pos(Laplas.parametric[i].svar,s);
            //ipos:=Find(s,Laplas.parametric[i].svar);

            // здесь корректно обрабатывается случай если две переменные
            // имеют похожие имена например $lg и $lgg. Если мы всречаем
            // $lgg то нужно подставить именно $lgg а не $lg. Так как вхождение обоих
            // начинается с одинаковой позиции.
            if ((bfirst or (ir=iposfirst)) and (ir>0)) then
            begin
               // если для одной и той же позиции в тексте
               // нашлось более одного вхождения то нужно
               // заменять вхождение наибольшей длины.
               bfirst:=false;
               iposfirst:=ir;
               ipos:=ir;
               if (ivar=-1) then
               begin
                  ivar:=i;
               end
               else
               begin
                  if (length(Laplas.parametric[i].svar)>length(Laplas.parametric[ivar].svar)) then
                  begin
                     ivar:=i;
                  end;
               end;
            end;
         end;
         if (ipos > 0) then
         begin
            // подстрока найдена, надо осуществить подстановку.
            bcont:=true; // успешная подстановка требуется продолжить поиск вхождений в дальнейшем.
            Delete(s,ipos,length(Laplas.parametric[ivar].svar));
            Insert('('+Laplas.parametric[ivar].sval+')',s,ipos);
         end;
      end;
   end;
   workstring:=s; // присваиваем преобразованную строку.
end;

procedure TFormVariables.all_obj_project_variable_rename(sold : String; snew : String);
 var
 i, i_4 : Integer;
 stext : String;
begin
    // блоки :
   for i:=0 to (Laplas.lb-1) do
   begin
      // цикл по всем блокам в расчётной области.


      // излучательная способность.
      stext:=Laplas.body[i].semissW;
      Laplas.body[i].semissW:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].semissE;
      Laplas.body[i].semissE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].semissS;
      Laplas.body[i].semissS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].semissN;
      Laplas.body[i].semissN:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].semissB;
      Laplas.body[i].semissB:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].semissT;
      Laplas.body[i].semissT:=StringReplace(stext,sold,snew,[rfReplaceAll]);

      stext:= Laplas.body[i].sxS;
      Laplas.body[i].sxS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.body[i].syS;
      Laplas.body[i].syS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.body[i].sxE;
      Laplas.body[i].sxE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.body[i].syE;
      Laplas.body[i].syE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.body[i].szS;
      Laplas.body[i].szS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.body[i].szE;
      Laplas.body[i].szE:=StringReplace(stext,sold,snew,[rfReplaceAll]);

      stext:=Laplas.body[i].sxC;
      Laplas.body[i].sxC:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].syC;
      Laplas.body[i].syC:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].szC;
      Laplas.body[i].szC:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].sHcyl;
      Laplas.body[i].sHcyl:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].sR_out_cyl;
      Laplas.body[i].sR_out_cyl:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:=Laplas.body[i].sR_in_cyl;
      Laplas.body[i].sR_in_cyl:=StringReplace(stext,sold,snew,[rfReplaceAll]);


     for i_4 := 0 to Laplas.body[i].n_power-1 do
     begin
        stext:= Laplas.body[i].arr_s_power[i_4];
        Laplas.body[i].arr_s_power[i_4]:=StringReplace(stext,sold,snew,[rfReplaceAll]);
     end;


   end;

   // источники тепла
   for i:=0 to (Laplas.ls-1) do
   begin
      // цикл по всем источникам тепла в расчётной области.

      // координаты плоского бесконечно тонкого источника тепла.
      stext:= Laplas.source[i].sxS;
      Laplas.source[i].sxS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.source[i].syS;
      Laplas.source[i].syS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.source[i].szS;
      Laplas.source[i].szS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.source[i].sxE;
      Laplas.source[i].sxE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.source[i].syE;
      Laplas.source[i].syE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.source[i].szE;
      Laplas.source[i].szE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.source[i].spower;
      Laplas.source[i].spower:=StringReplace(stext,sold,snew,[rfReplaceAll]);

   end;

   // твёрдые стенки
   for i:=0 to (Laplas.lw-1) do
   begin
      // цикл по всем источникам тепла в расчётной области.


      // координаты стенки
      stext:= Laplas.wall[i].sxS;
      Laplas.wall[i].sxS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.wall[i].syS;
      Laplas.wall[i].syS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.wall[i].szS;
      Laplas.wall[i].szS:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.wall[i].sxE;
      Laplas.wall[i].sxE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.wall[i].syE;
      Laplas.wall[i].syE:=StringReplace(stext,sold,snew,[rfReplaceAll]);
      stext:= Laplas.wall[i].szE;
      Laplas.wall[i].szE:=StringReplace(stext,sold,snew,[rfReplaceAll]);

       stext:= Laplas.wall[i].semissivity;
       Laplas.wall[i].semissivity:=StringReplace(stext,sold,snew,[rfReplaceAll]);
       stext:= Laplas.wall[i].sheat_transfer_coefficient;
       Laplas.wall[i].sheat_transfer_coefficient:=StringReplace(stext,sold,snew,[rfReplaceAll]);

      end;
end;

// рекурсивная функция конструирования поддерева
// выражения с корнем r из строки st
procedure TFormVariables.ButtonRenameClick(Sender: TObject);
var
 i : Integer;
 sold, snew : String;
begin
    // Показать форму в которой надо выбрать что мы удаляем.
    FormRenameVar.ComboBox1.Clear();
    for i := 0 to Laplas.ivar-1 do
      begin
         FormRenameVar.ComboBox1.Items.Add(Trim(Laplas.parametric[i].svar));
      end;
      FormRenameVar.ComboBox1.ItemIndex:=Laplas.ivar-1;
      FormRenameVar.EditNewName.Text:='';
    FormRenameVar.ShowModal;
    if (FormRenameVar.bOk_rename) then
    begin
       snew:= Trim(FormRenameVar.EditNewName.Text);
       sold:=Trim(FormRenameVar.ComboBox1.Items[FormRenameVar.ComboBox1.ItemIndex]);

       all_obj_project_variable_rename(sold,snew);

       FormVariables.StringGridVariables.Cells[1,1+FormRenameVar.ComboBox1.ItemIndex]:=Trim(snew);
       Laplas.parametric[FormRenameVar.ComboBox1.ItemIndex].svar:=Trim(snew);

       // обновление размеров объектов
       my_update_size();
    end;
end;

procedure TFormVariables.Constr_Tree(var r : pTop; var st : string);
var
    po,code : Integer;
    stl, stri : String;
    c : Real;
    bOk : Boolean;
    sbuf : String; // временная строка для корректного введения новой переменной.
    bcont : Boolean;
   stxe8 : String;
   rootbuf : pTop;

// внутренняя функция поиска разделительного знака в строке st:
// SetOp - множество знаков; функция возвращает позицию разделительного знака или 0.
function PosOp(st : String; SetOp : setChar) : byte;
var i,j,k,p : byte;
begin
   j:=0; k:=0; p:=0;
   i:=length(st);//i:=1; // Внимание !!! строку обязательно нужно сканировать из конца в начало.
   while ((i>=1) and (i<=length(st)) and (p=0)) do
   begin
      if st[i]='(' then inc(j) // считаем количество открывающихся скобок
      else if st[i]=')' then inc(k) // считаем количество закрывающихся скобок
      else if (j=k) and (st[i] in SetOp) then p:=i;
      dec(i); //inc(i); // Внимание строку обязательно нужно сканировать из конца в начало.
   end;
   PosOp:=p;
end;

// раздел операторов функции конструирования дерева выражения

begin
   po:=PosOp(st,['+','-']); // ищем разделительный знак оперпации + или -
   if (po<>0) then
   begin
      val(StringReplace(Trim(st),',','.',[rfReplaceAll]),c,code);
      if code=0 then
      begin
         po:=0; // случай 75.0e-3
      end;
   end;
   if po=0 then po:=PosOp(st,['*','/']); // ищем разделительный знак операции  * или /
   if po=0 then po:=PosOp(st,['^']); // ищем разделительный знак операции возведения в степень
   if (po<>0) then  // разделяющий знак найден
   begin
      stxe8:=' ';
      stxe8[1]:=st[po];
      //(r^).myoperator[1]:=st[po]; // записываем знак операции в вершину
      (r^).myoperator:=stxe8;

      stl:=copy(st,1,po-1); // копируем подстроку первого операнда
      if (length(stl)>0) then
      begin
         if (stl[1]='(') and (PosOp(stl,['*','/','+','-','^'])=0) then
                      stl:=copy(stl,2,length(stl)-2); // убираем скобки
      end
       else
      begin
         // унарный минус
         stl:='';
         bOk:=false;
         if (((r^).myoperator[1]='-') or (r^.myoperator[1]='+')) then
         begin
            bOk:=true;  //b1:=true;
         end;
         if (not bOk) then
         begin
            // дерево конечно построится но вычислить значение выражения по
            // нему будет нельзя. Этот случай нужно корректно обработать, т.е.
            // сообщить об ошибке наверх по иерархии вызовывов, что бы избежать
            // вызова приводящего к ошибке вычисления выражения по такому днреву.
            Laplas.MainMemo.Lines.Add('Error! tree can not be construct... ');
            Application.MessageBox(PChar('tree can not be construct...'),'Error!',MB_OK);
         end;
      end;

      stri:=copy(st,po+1,length(st)-po); // копируем подстроку второго операнда
      if (stri[1]='(') and (PosOp(stri,['*','/','+','-','^'])=0) then
      begin
          stri:=copy(stri,2,length(stri)-2); // убираем скобки
          //b1:=false; // Особый случай снят с рассмотрения. У нас был минус перед большой скобкой.
      end;


      if (length(stl)>0) then
      begin
         rootbuf:=nil;
         new(rootbuf);
         //new(r^.left);  // создаём левое поддерево
         //Constr_Tree(r^.left,stl); // конструируем левый операнд
         Constr_Tree(rootbuf,stl);
         r^.left:=rootbuf;
         rootbuf:=nil;

      end
       else
      begin
         // унарный минус
         r^.left:=nil; // левое поддерево пустое.
      end;
      rootbuf:=nil;
      new(rootbuf);
      //new(r^.right); // создаём правое поддерево
      //Constr_Tree(r^.right,stri); // конструируем правый операнд
      Constr_Tree(rootbuf,stri);
      r^.right:=rootbuf;
      rootbuf:=nil;


    end
    else
    begin
       // Здесь предполагается что все переменные были заранее исключены
       // применением процедуры  my_substitutional_of_value.
       bcont:=True;
       while (bcont) do
       begin
          st:=Trim(st); // Убираем пробелы.
          if (length(st)>0) then
          begin
             if ((st[1]='(') and (st[Length(st)]=')')) then
             begin
                st:=copy(st,2,length(st)-2); // убираем скобки
             end
              else
             begin
                bcont:=False;
             end;
          end;
       end;

       if (length(st)>0) then
       begin
           val(StringReplace(st,',','.',[rfReplaceAll]),c,code); // пытаемся получить число
       end
       else
       begin
          // Пустая строка.
          code:=0;
          c:=0.0;
       end;

       if code=0 then // константа
       begin
          r^.myoperator:='o';
          r^.left:=nil;
          r^.right:=nil;
          r^.value:=c;
       end
       else // функция
       begin
          po:=Pos('(',st);
          if (po=0) then
          begin
             sbuf:=st;
             my_substitutional_of_value(sbuf,st); // результат записан в st.
             // Теперь если переменная уже была введена, то st будет просто
             //  константой и вводить новую переменную не нужно.
             val(StringReplace(st,',','.',[rfReplaceAll]),c,code); // пытаемся получить число
             if (code=0) then
             begin
                // Ранее по ходу построения дерева эта переменная уже была введена,
                // и поэтому теперь мы записываем просто её константное значение.

                r^.myoperator:='o';
                r^.left:=nil;
                r^.right:=nil;
                r^.value:=c;
             end
              else
             begin
                // Это точно новая переменная.
                AddVariableForm.lblname.Caption:=Trim(st); // предполагаемое имя переменной.
                if (FormatSettings.DecimalSeparator='.') then
                begin
                   AddVariableForm.edtvalue.Text:='0.0'; // предполагаемое значение.
                end;
                if (FormatSettings.DecimalSeparator=',') then
                begin
                   AddVariableForm.edtvalue.Text:='0,0'; // предполагаемое значение.
                end;
                AddVariableForm.ShowModal;
                // переменная воспринимается как константное значение этой переменной.
                r^.myoperator:='o';
                r^.left:=nil;
                r^.right:=nil;

                val(StringReplace(AddVariableForm.edtvalue.Text,',','.',[rfReplaceAll]),c,code);
                if (code=0) then
                begin
                   r^.value:=c;
                end
                 else
                begin
                   // Данного случая просто не может быть, см модуль AddVariableUnit !!!
                end;
             end;
          end
           else
          begin
             r^.myoperator:=copy(st,1,po-1); // выделяем имя функции
             r^.right:=nil;
             stl:=copy(st,po+1,length(st)-po-1); // выделяем подстроку параметра
             rootbuf:=nil;
             new(rootbuf);

             //new(r^.left);
             //Constr_Tree(r^.left,stl); // конструируем параметр
             Constr_Tree(rootbuf,stl); // конструируем параметр
             r^.left:=rootbuf;
             rootbuf:=nil;
          end;
       end;
    end;
end;   // Constr_Tree

// рекурсивное вычисление значения функции
// если key=false, то значение не существует
function TFormVariables.my_count(r : pTop; var key : Boolean) : Real;
var
    s,s1 : Real;
begin
   if not key then // значение функции не существует
   begin
      my_count:=0.0;
      Laplas.MainMemo.Lines.Add('Error! expression can not be calculated... ');
      Application.MessageBox(PChar('expression can not be calculated...'),'Error!',MB_OK);
   end
   else
   begin
      if r^.myoperator[1]='o' then my_count:=r^.value // константа
      else
       case r^.myoperator[1] of
         '+' :  begin
                    if (r^.left=nil) then
                   begin
                      // унарный плюс (есть вероятность что он может встретиться.
                      my_count:=my_count(r^.right,key);
                   end
                   else
                   begin
                      my_count:=my_count(r^.left,key)+my_count(r^.right,key);
                   end;
                end;
         '-' :  begin
                   if (r^.left=nil) then
                   begin
                      // унарный минус
                      my_count:=-my_count(r^.right,key);
                   end
                   else
                   begin
                      my_count:=my_count(r^.left,key)-my_count(r^.right,key);
                   end;
                end;
         '*' :  my_count:=my_count(r^.left,key)*my_count(r^.right,key);
         '/' :  begin
                   s:=my_count(r^.right,key);
                   if abs(s)<1e-10 then // практический ноль
                   begin
                       my_count:=0.0;
                       key:=false;
                   end
                   else my_count:=my_count(r^.left,key)/s;
                end;
         '^' : begin
                  s:=my_count(r^.left,key);
                  s1:=my_count(r^.right,key);
                  if abs(s)<1e-10 then // практический ноль
                  begin
                     if abs(s1)<1e-10 then my_count:=1.0
                     else  my_count:=0.0;
                  end
                  else
                  begin
                      my_count:=exp(s1*ln(abs(s)));
                  end;
               end;
         // некоторые немногие основные математические функции
         'a' : begin
                  // abs
                  my_count:=abs(my_count(r^.left,key));
               end;
         's' : begin
                  if  r^.myoperator[2]='i' then
                  begin
                     // sin
                     my_count:=sin(my_count(r^.left,key));
                  end;
                  if  r^.myoperator[2]='h' then
                  begin
                     // sh - гиперболический синус
                     s:=my_count(r^.left,key);
                     my_count:=0.5*(exp(s)-exp(-s));
                  end;
                  if  r^.myoperator[2]='q' then
                  begin
                     if (length(r^.myoperator[2])=4) then
                     begin
                        // sqrt(abs(s)) квадратный корень
                        my_count:=sqrt(abs(my_count(r^.left,key)));
                     end
                     else
                     begin
                        // sqr - возведение в квадрат
                        my_count:=sqr(my_count(r^.left,key));
                     end;
                  end;
               end;
         'c' : begin
                  if  r^.myoperator[2]='o' then
                  begin
                     // cos
                     my_count:=cos(my_count(r^.left,key));
                  end;
                  if  r^.myoperator[2]='h' then
                  begin
                     // ch гиперболический косинус
                     s:=my_count(r^.left,key);
                     my_count:=0.5*(exp(s)+exp(-s));
                  end;
               end;
         'e' : begin
                   // exp   экспонента
                   my_count:=exp(my_count(r^.left,key));
               end;
         'l' : begin
                  // ln(abs(s))  натуральный логарифм
                  s:=my_count(r^.left,key);
                  my_count:=ln(abs(s));
               end;
         else // неопределённая операция
          begin
             my_count:=0.0;
             key:=false;
          end;
       end; // case

   end;
end; // my_count

// освобождение оперативной памяти занимаемой двоичным деревом
procedure TFormVariables.my_delete(var r : pTop);
var
 rootbuf : pTop;
begin
    if (r<>nil) then
    begin
        rootbuf:=nil;
        rootbuf:=r^.left;
        r^.left:=nil;
        //my_delete(r^.left);
        my_delete(rootbuf);
        rootbuf:=nil;
        rootbuf:=r^.right;
        r^.right:=nil;
        //my_delete(r^.right);
        my_delete(rootbuf);
        rootbuf:=nil;
        Dispose(r);
    end;
end;  // my_delete

// преобразует параметризованую строку содержащую переменные в вещественное число.
function TFormVariables.my_real_convert(s : String; var bOk : Boolean) : Real;
var
   s1 : String;  // рабочая строка
   Root : pTop; // корень дерева выражения
   key : Boolean; // признак существования выражения в заданной точке
   r1 : Real; // возвращаемое значение
   
begin
    // подстановка - замена всех переменных их значениями:
    my_substitutional_of_value(s ,s1);
    // Если переменные остались, то их необходимо добавить в список пользовательских переменных.
    // построение дерева.
    new(Root); // выделение оперативной памяти
    Constr_Tree(Root, s1); // построение дерева выражения
    key:=true;
    r1:=my_count(Root, key);
    my_delete(Root); // освобождение оперативной памяти
    bOk:=key; // передаваемый наверх признак успешности операции.
    if not key then
    begin
       r1:=0.0;
       Laplas.MainMemo.Lines.Add('Error! expression :');
       Laplas.MainMemo.Lines.Add(s);
       Laplas.MainMemo.Lines.Add('can not be calculated... ');
       Application.MessageBox(PChar('expression '+s+' can not be calculated...'),'Error!',MB_OK);
    end;
    my_real_convert:=r1;
end;

// обновление размеров после изменения значения переменной
// которая задана пользователем.
procedure TFormVariables.my_update_size();
var
    i, i_4 : Integer;
    // вспомогательные переменные для обработки исключительной ситуации
    bOk : Boolean;
    r1, r2, r3, r4, r5, r6, rpow : Real;
    r7, r8, r9, r10, r11, r12 : Real;

begin
   // инициализация :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;
   r4:=0.0;
   r5:=0.0;
   r6:=0.0;
   rpow:=0.0;
   r7:=0.0;
   r8:=0.0;
   r9:=0.0;
   r10:=0.0;
   r11:=0.0;
   r12:=0.0;

   // блоки :
   for i:=0 to (Laplas.lb-1) do
   begin
      // цикл по всем блокам в расчётной области.
      bOk:=true; // признак правильности ввода

      // излучательная способность.

      if bOk then r1:=FormVariables.my_real_convert(Laplas.body[i].semissW,bOk);  // числовые размеры
      if bOk then r2:=FormVariables.my_real_convert(Laplas.body[i].semissE,bOk);  // заданные пользователем
      if bOk then r3:=FormVariables.my_real_convert(Laplas.body[i].semissS,bOk);  // с учётом подстановки
      if bOk then r4:=FormVariables.my_real_convert(Laplas.body[i].semissN,bOk);  // значений переменных.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.body[i].semissB,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.body[i].semissT,bOk);

      if (bOk)  then
      begin
        // Обновление значений излучательной способности.
       Laplas.body[i].emissW:=r1;
       Laplas.body[i].emissE:=r2;
       Laplas.body[i].emissS:=r3;
       Laplas.body[i].emissN:=r4;
       Laplas.body[i].emissB:=r5;
       Laplas.body[i].emissT:=r6;
      end;


      // координаты блока

      if bOk then r1:=FormVariables.my_real_convert(Laplas.body[i].sxS,bOk);  // числовые размеры
      if bOk then r2:=FormVariables.my_real_convert(Laplas.body[i].syS,bOk);  // заданные пользователем
      if bOk then r3:=FormVariables.my_real_convert(Laplas.body[i].sxE,bOk);  // с учётом подстановки
      if bOk then r4:=FormVariables.my_real_convert(Laplas.body[i].syE,bOk);  // значений переменных.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.body[i].szS,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.body[i].szE,bOk);

        if bOk then r7:=FormVariables.my_real_convert(Laplas.body[i].sxC,bOk);  // числовые размеры
      if bOk then r8:=FormVariables.my_real_convert(Laplas.body[i].syC,bOk);  // заданные пользователем
      if bOk then r9:=FormVariables.my_real_convert(Laplas.body[i].szC,bOk);  // с учётом подстановки
      if bOk then r10:=FormVariables.my_real_convert(Laplas.body[i].sHcyl,bOk);  // значений переменных.
      if bOk then r11:=FormVariables.my_real_convert(Laplas.body[i].sR_out_cyl,bOk);
      if bOk then r12:=FormVariables.my_real_convert(Laplas.body[i].sR_in_cyl,bOk);
      if bOk then
      begin
         for i_4 := 0 to Laplas.body[i].n_power-1 do
         begin
            if (bOk) then
            begin
               rpow:=FormVariables.my_real_convert(Laplas.body[i].arr_s_power[i_4],bOk);
            end;
         end;
      end;

      if (bOk) then
      begin
         Laplas.body[i].xS:=r1;  // числовые размеры
         Laplas.body[i].yS:=r2;  // заданные пользователем
         Laplas.body[i].xE:=r3;  // с учётом подстановки
         Laplas.body[i].yE:=r4;  // значений переменных.
         Laplas.body[i].zS:=r5;
         Laplas.body[i].zE:=r6;

         Laplas.body[i].xC:=r7;  // числовые размеры
         Laplas.body[i].yC:=r8;  // заданные пользователем
         Laplas.body[i].zC:=r9;  // с учётом подстановки
         Laplas.body[i].Hcyl:=r10;  // значений переменных.
         Laplas.body[i].R_out_cyl:=r11;
         Laplas.body[i].R_in_cyl:=r12;
         for i_4 := 0 to Laplas.body[i].n_power-1 do
         begin
             Laplas.body[i].arr_power[i_4]:=FormVariables.my_real_convert(Laplas.body[i].arr_s_power[i_4],bOk); // мощность тепловыделения
         end;
      end;

   end;

   // источники тепла
   for i:=0 to (Laplas.ls-1) do
   begin
      // цикл по всем источникам тепла в расчётной области.
      bOk:=true; // признак правильности ввода

      // координаты блока

      if bOk then r1:=FormVariables.my_real_convert(Laplas.source[i].sxS,bOk);  // числовые размеры
      if bOk then r2:=FormVariables.my_real_convert(Laplas.source[i].syS,bOk);  // заданные пользователем
      if bOk then r3:=FormVariables.my_real_convert(Laplas.source[i].sxE,bOk);  // с учётом подстановки
      if bOk then r4:=FormVariables.my_real_convert(Laplas.source[i].syE,bOk);  // значений переменных.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.source[i].szS,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.source[i].szE,bOk);
      if bOk then rpow:=FormVariables.my_real_convert(Laplas.source[i].spower,bOk);

      if (bOk) then
      begin
         Laplas.source[i].xS:=r1;  // числовые размеры
         Laplas.source[i].yS:=r2;  // заданные пользователем
         Laplas.source[i].xE:=r3;  // с учётом подстановки
         Laplas.source[i].yE:=r4;  // значений переменных.
         Laplas.source[i].zS:=r5;
         Laplas.source[i].zE:=r6;
         Laplas.source[i].Power:=rpow; // мощность тепловыделения
      end;
   end;

   // твёрдые стенки
   for i:=0 to (Laplas.lw-1) do
   begin
      // цикл по всем источникам тепла в расчётной области.
      bOk:=true; // признак правильности ввода

      // координаты блока

      if bOk then r1:=FormVariables.my_real_convert(Laplas.wall[i].sxS,bOk);  // числовые размеры
      if bOk then r2:=FormVariables.my_real_convert(Laplas.wall[i].syS,bOk);  // заданные пользователем
      if bOk then r3:=FormVariables.my_real_convert(Laplas.wall[i].sxE,bOk);  // с учётом подстановки
      if bOk then r4:=FormVariables.my_real_convert(Laplas.wall[i].syE,bOk);  // значений переменных.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.wall[i].szS,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.wall[i].szE,bOk);

      if (bOk) then
      begin
         Laplas.wall[i].xS:=r1;  // числовые размеры
         Laplas.wall[i].yS:=r2;  // заданные пользователем
         Laplas.wall[i].xE:=r3;  // с учётом подстановки
         Laplas.wall[i].yE:=r4;  // значений переменных.
         Laplas.wall[i].zS:=r5;
         Laplas.wall[i].zE:=r6;
      end;

       bOk:=true; // признак правильности ввода
       if bOk then r1:=FormVariables.my_real_convert(Laplas.wall[i].semissivity,bOk);  // числовые размеры
       if bOk then r2:=FormVariables.my_real_convert(Laplas.wall[i].sheat_transfer_coefficient,bOk);

      if (bOk) then
      begin
         // излучательная способность.
         Laplas.wall[i].emissivity:=r1;
         // коэффициент теплоотдачи.
         Laplas.wall[i].heat_transfer_coefficient:=r2;
      end;

   end;

end;

end.
