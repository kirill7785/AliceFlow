unit UnitVariables;
// ���������� ��� �������� ������������ � �������� Optimetric`a.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, Grids, StdCtrls;

type

  // ���� ��� ���������� ������������������ ���������
  setChar=set of char; // ��� ��������� ��������
  pTop=^Top;
  Top=record
  private
     //myoperator : string[5];
     Fmyoperator : string; // ���� ��������
     Fvalue : Real; // �������� ���������
     Fleft, Fright : pTop; // ��������� �� ����� � ������ ���������
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

      // ����� ��������� � ������
      //function Find(const S, P : String) : Integer;
      // ����������� �������� ���������� � ������
      // patterns - �������� ������������ ������.
      // workstring - ��������������� �������� ������.
      procedure my_substitutional_of_value(patterns : String; var workstring : String);
      // ����������� ������� ��������������� ���������
      // ��������� � ������ r �� ������ st
      procedure Constr_Tree(var r : pTop; var st : string);
      // ����������� ���������� �������� �������
      // ���� key=false, �� �������� �� ����������
      function my_count(r : pTop;  var key : Boolean) : Real;
      // ������������ ����������� ������ ���������� �������� �������
      procedure my_delete(var r : pTop);



  public
    { Public declarations }

      // ����������� ���������������� ������ ���������� ���������� � ������������ �����.
      // ���� bOk = true �� �������� �������������� ������ �������, ���� bOk=false
      // �� �������� �������������� �������� ������ � �������� �������������� ���� ��������
      function my_real_convert(s : String; var bOk : Boolean) : Real;
      procedure all_obj_project_variable_rename(sold : String; snew : String);
      // ���������� �������� ����� ��������� �������� ����������
      // ������� ������ �������������.
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

// ������������� ����� ����� ����� ��������.
procedure TFormVariables.FormCreate(Sender: TObject);
var
    i : Integer;
begin
   // ������������� ����� �������������� ��� ����������� :
   // Cells[�������, ������]
   FormVariables.StringGridVariables.Cells[1,0]:='$var';
   FormVariables.StringGridVariables.Cells[2,0]:='value';
   for i:=1 to FormVariables.StringGridVariables.RowCount-1 do
   begin
      FormVariables.StringGridVariables.Cells[0,i]:=IntToStr(i);
   end;

    for i:=1 to FormVariables.StringGridVariables.RowCount-1 do
    begin
       if (i-1>=Laplas.ivar) then
       begin
          FormVariables.StringGridVariables.Cells[1,i]:='';
          FormVariables.StringGridVariables.Cells[2,i]:='';
       end
       else
       begin
          FormVariables.StringGridVariables.Cells[1,i]:=Laplas.parametric[i-1].svar;
          FormVariables.StringGridVariables.Cells[2,i]:=Laplas.parametric[i-1].sval;
       end;
    end;
end;

procedure TFormVariables.BApplyClick(Sender: TObject);
var
     i : Integer;
     b : Boolean;
     s : String;
begin
    // ������ ���������� � ������.
    // ���������� ������ ���� ��������������� ������ ����.
    // ���������� �� ��������������� ����� �������� ����������� �� ������������.
    // ��� ���������� ����������� ������ ���������� �� ����� $.
    b:=true;
    i:=0;
    while (b and (i<StringGridVariables.RowCount-2)) do
    begin
       s:=StringGridVariables.Cells[1,i+1];
       if ((length(s)>=2) and (s[1]='$')) then b:=true
       else b:=false;
       if (b) then i:=i+1; // ��������� � ��������� ����������
    end;
    Laplas.ivar:=i; // ���������� ����������.
    SetLength(Laplas.parametric,Laplas.ivar);
    for i:=0 to Laplas.ivar-1 do
    begin
       Laplas.parametric[i].svar:=Trim(StringGridVariables.Cells[1,i+1]);
       Laplas.parametric[i].sval:=Trim(StringGridVariables.Cells[2,i+1]);
    end;
    // ��� ���������� �������� � �������  StringGridVariables

    // ��� ������������ ���������� ������������ � �� ���� �������
    for i:=Laplas.ivar+1 to StringGridVariables.RowCount-1 do
    begin
       StringGridVariables.Cells[1,i]:='';
       StringGridVariables.Cells[2,i]:='';
    end;

    // ���������� �������� ��������
    my_update_size();
    // ���������� ���������
    Laplas.ReadyPaint;
    Close();
end;

// ����� ��������� � ������
// ������ ���������� ������� pos.
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

// ����������� �������� ���������� � ������
// patterns - �������� ������������ ������.
// workstring - ��������������� �������� ������.
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
         // �������������
         bcont:=false;
         ipos:=0;
         bfirst:=true;
         iposfirst:=0;
         ivar:=-1; // �������������� ����������
         for i:=0 to Laplas.ivar-1 do
         begin
            ir:=Pos(Laplas.parametric[i].svar,s);
            //ipos:=Find(s,Laplas.parametric[i].svar);

            // ����� ��������� �������������� ������ ���� ��� ����������
            // ����� ������� ����� �������� $lg � $lgg. ���� �� ��������
            // $lgg �� ����� ���������� ������ $lgg � �� $lg. ��� ��� ��������� �����
            // ���������� � ���������� �������.
            if ((bfirst or (ir=iposfirst)) and (ir>0)) then
            begin
               // ���� ��� ����� � ��� �� ������� � ������
               // ������� ����� ������ ��������� �� �����
               // �������� ��������� ���������� �����.
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
            // ��������� �������, ���� ����������� �����������.
            bcont:=true; // �������� ����������� ��������� ���������� ����� ��������� � ����������.
            Delete(s,ipos,length(Laplas.parametric[ivar].svar));
            Insert('('+Laplas.parametric[ivar].sval+')',s,ipos);
         end;
      end;
   end;
   workstring:=s; // ����������� ��������������� ������.
end;

procedure TFormVariables.all_obj_project_variable_rename(sold : String; snew : String);
 var
 i, i_4 : Integer;
 stext : String;
begin
    // ����� :
   for i:=0 to (Laplas.lb-1) do
   begin
      // ���� �� ���� ������ � ��������� �������.


      // ������������� �����������.
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

   // ��������� �����
   for i:=0 to (Laplas.ls-1) do
   begin
      // ���� �� ���� ���������� ����� � ��������� �������.

      // ���������� �������� ���������� ������� ��������� �����.
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

   // ������ ������
   for i:=0 to (Laplas.lw-1) do
   begin
      // ���� �� ���� ���������� ����� � ��������� �������.


      // ���������� ������
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

// ����������� ������� ��������������� ���������
// ��������� � ������ r �� ������ st
procedure TFormVariables.ButtonRenameClick(Sender: TObject);
var
 i : Integer;
 sold, snew : String;
begin
    // �������� ����� � ������� ���� ������� ��� �� �������.
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

       // ���������� �������� ��������
       my_update_size();
    end;
end;

procedure TFormVariables.Constr_Tree(var r : pTop; var st : string);
var
    po,code : Integer;
    stl, stri : String;
    c : Real;
    bOk : Boolean;
    sbuf : String; // ��������� ������ ��� ����������� �������� ����� ����������.
    bcont : Boolean;
   stxe8 : String;
   rootbuf : pTop;

// ���������� ������� ������ ��������������� ����� � ������ st:
// SetOp - ��������� ������; ������� ���������� ������� ��������������� ����� ��� 0.
function PosOp(st : String; SetOp : setChar) : byte;
var i,j,k,p : byte;
begin
   j:=0; k:=0; p:=0;
   i:=length(st);//i:=1; // �������� !!! ������ ����������� ����� ����������� �� ����� � ������.
   while ((i>=1) and (i<=length(st)) and (p=0)) do
   begin
      if st[i]='(' then inc(j) // ������� ���������� ������������� ������
      else if st[i]=')' then inc(k) // ������� ���������� ������������� ������
      else if (j=k) and (st[i] in SetOp) then p:=i;
      dec(i); //inc(i); // �������� ������ ����������� ����� ����������� �� ����� � ������.
   end;
   PosOp:=p;
end;

// ������ ���������� ������� ��������������� ������ ���������

begin
   st:=Trim(st);
   //Application.MessageBox(PWideChar(st),'err',MB_OK);
   while (st[1]='(') and (PosOp(st,['*','/','+','-','^'])=0) do
   begin
                      st:=copy(st,2,length(st)-2);
                      st:=Trim(st);
                     // Application.MessageBox(PWideChar(st),'err',MB_OK);
   end;
   po:=PosOp(st,['+','-']); // ���� �������������� ���� ��������� + ��� -
   if (po<>0) then
   begin
      if ((po>2)and((st[po-1]='e')or(st[po-1]='E'))
      and(((st[po-2]>='0')and(st[po-2]<='9'))or(st[po-2]=',')or(st[po-2]='.'))
      and(po+1<=length(st))
      and(st[po+1]>='0')and(st[po+1]<='9')) then
      begin
         // ���������������� ������ �����.
         po:=0;
      end
      else
      begin
         if (FormatSettings.DecimalSeparator='.') then
         begin
            // ����������� ����� � ������� ����� �����.
            // ������ ��� ������� �� �����.
            val(StringReplace(Trim(st),',','.',[rfReplaceAll]),c,code);
         end
         else
         begin
            // ����������� ����� � ������� ����� �������.
            // ������ ��� ����� �� �������.
            val(StringReplace(Trim(st),'.',',',[rfReplaceAll]),c,code);
         end;
         if code=0 then
         begin
            po:=0; // ������ 75.0e-3
         end;
      end;
   end;
   if po=0 then po:=PosOp(st,['*','/']); // ���� �������������� ���� ��������  * ��� /
   if po=0 then po:=PosOp(st,['^']); // ���� �������������� ���� �������� ���������� � �������
   if (po<>0) then  // ����������� ���� ������
   begin
      stxe8:=' ';
      stxe8[1]:=st[po];
      //(r^).myoperator[1]:=st[po]; // ���������� ���� �������� � �������
      (r^).myoperator:=stxe8;

      stl:=Trim(copy(st,1,po-1)); // �������� ��������� ������� ��������
      if (length(stl)>0) then
      begin
         //Application.MessageBox(PWideChar(stl),'err',MB_OK);
         while (stl[1]='(') and (PosOp(stl,['*','/','+','-','^'])=0) do
         begin
            stl:=Trim(copy(stl,2,length(stl)-2)); // ������� ������

         //Application.MessageBox(PWideChar(stl),'err',MB_OK);
         end;
      end
       else
      begin
         // ������� �����
         stl:='';
         bOk:=false;
         if (((r^).myoperator[1]='-') or (r^.myoperator[1]='+')) then
         begin
            bOk:=true;  //b1:=true;
         end;
         if (not bOk) then
         begin
            // ������ ������� ���������� �� ��������� �������� ��������� ��
            // ���� ����� ������. ���� ������ ����� ��������� ����������, �.�.
            // �������� �� ������ ������ �� �������� ���������, ��� �� ��������
            // ������ ����������� � ������ ���������� ��������� �� ������ ������.
            Laplas.MainMemo.Lines.Add('Error! tree can not be construct... ');
            Application.MessageBox(PChar('tree can not be construct...'),'Error!',MB_OK);
         end;
      end;


      stri:=Trim(copy(st,po+1,length(st)-po)); // �������� ��������� ������� ��������
      //Application.MessageBox(PWideChar(stri),'err',MB_OK);
      while (stri[1]='(') and (PosOp(stri,['*','/','+','-','^'])=0) do
      begin
          stri:=Trim(copy(stri,2,length(stri)-2)); // ������� ������
          //b1:=false; // ������ ������ ���� � ������������. � ��� ��� ����� ����� ������� �������.
          //Application.MessageBox(PWideChar(stri),'err',MB_OK);
      end;


      if (length(stl)>0) then
      begin
         rootbuf:=nil;
         new(rootbuf);
         //new(r^.left);  // ������ ����� ���������
         //Constr_Tree(r^.left,stl); // ������������ ����� �������
         Constr_Tree(rootbuf,stl);
         r^.left:=rootbuf;
         rootbuf:=nil;

      end
       else
      begin
         // ������� �����
         r^.left:=nil; // ����� ��������� ������.
      end;
      rootbuf:=nil;
      new(rootbuf);
      //new(r^.right); // ������ ������ ���������
      //Constr_Tree(r^.right,stri); // ������������ ������ �������
      Constr_Tree(rootbuf,stri);
      r^.right:=rootbuf;
      rootbuf:=nil;


    end
    else
    begin
       // ����� �������������� ��� ��� ���������� ���� ������� ���������
       // ����������� ���������  my_substitutional_of_value.
       bcont:=True;
       while (bcont) do
       begin
          st:=Trim(st); // ������� �������.
          if (length(st)>0) then
          begin
             if ((st[1]='(') and (st[Length(st)]=')')) then
             begin
                st:=copy(st,2,length(st)-2); // ������� ������
             end
              else
             begin
                bcont:=False;
             end;
          end;
       end;

       if (length(st)>0) then
       begin
       if (FormatSettings.DecimalSeparator='.') then
         begin
            // ����������� ����� � ������� ����� �����.
            // ������ ��� ������� �� �����.
            val(StringReplace(Trim(st),',','.',[rfReplaceAll]),c,code);  // �������� �������� �����
         end
         else
         begin
            // ����������� ����� � ������� ����� �������.
            // ������ ��� ����� �� �������.
            val(StringReplace(Trim(st),'.',',',[rfReplaceAll]),c,code);  // �������� �������� �����
         end;
       end
       else
       begin
          // ������ ������.
          code:=0;
          c:=0.0;
       end;

       if code=0 then // ���������
       begin
          r^.myoperator:='o';
          r^.left:=nil;
          r^.right:=nil;
          r^.value:=c;
       end
       else // �������
       begin
          po:=Pos('(',st);
          if (po=0) then
          begin
             sbuf:=st;
             my_substitutional_of_value(sbuf,st); // ��������� ������� � st.
             // ������ ���� ���������� ��� ���� �������, �� st ����� ������
             //  ���������� � ������� ����� ���������� �� �����.
             if (FormatSettings.DecimalSeparator='.') then
             begin
                // ����������� ����� � ������� ����� �����.
                // ������ ��� ������� �� �����.
                val(StringReplace(Trim(st),',','.',[rfReplaceAll]),c,code);  // �������� �������� �����
             end
              else
             begin
                // ����������� ����� � ������� ����� �������.
                // ������ ��� ����� �� �������.
                val(StringReplace(Trim(st),'.',',',[rfReplaceAll]),c,code);  // �������� �������� �����
             end;

             if (code=0) then
             begin
                // ����� �� ���� ���������� ������ ��� ���������� ��� ���� �������,
                // � ������� ������ �� ���������� ������ � ����������� ��������.

                r^.myoperator:='o';
                r^.left:=nil;
                r^.right:=nil;
                r^.value:=c;
             end
              else
             begin
                // ��� ����� ����� ����������.
                AddVariableForm.lblname.Caption:=Trim(st); // �������������� ��� ����������.
                if (FormatSettings.DecimalSeparator='.') then
                begin
                   AddVariableForm.edtvalue.Text:='0.0'; // �������������� ��������.
                end;
                if (FormatSettings.DecimalSeparator=',') then
                begin
                   AddVariableForm.edtvalue.Text:='0,0'; // �������������� ��������.
                end;
                AddVariableForm.ShowModal;
                // ���������� �������������� ��� ����������� �������� ���� ����������.
                r^.myoperator:='o';
                r^.left:=nil;
                r^.right:=nil;

                if (FormatSettings.DecimalSeparator='.') then
                begin
                   // ����������� ����� � ������� ����� �����.
                   // ������ ��� ������� �� �����.
                   val(StringReplace(AddVariableForm.edtvalue.Text,',','.',[rfReplaceAll]),c,code);  // �������� �������� �����
                end
                 else
                begin
                   // ����������� ����� � ������� ����� �������.
                   // ������ ��� ����� �� �������.
                   val(StringReplace(AddVariableForm.edtvalue.Text,'.',',',[rfReplaceAll]),c,code);  // �������� �������� �����
                end;

                if (code=0) then
                begin
                   r^.value:=c;
                end
                 else
                begin
                   // ������� ������ ������ �� ����� ����, �� ������ AddVariableUnit !!!
                end;
             end;
          end
           else
          begin
             r^.myoperator:=copy(st,1,po-1); // �������� ��� �������
             r^.right:=nil;
             stl:=copy(st,po+1,length(st)-po-1); // �������� ��������� ���������
             rootbuf:=nil;
             new(rootbuf);

             //new(r^.left);
             //Constr_Tree(r^.left,stl); // ������������ ��������
             Constr_Tree(rootbuf,stl); // ������������ ��������
             r^.left:=rootbuf;
             rootbuf:=nil;
          end;
       end;
    end;
end;   // Constr_Tree

// ����������� ���������� �������� �������
// ���� key=false, �� �������� �� ����������
function TFormVariables.my_count(r : pTop; var key : Boolean) : Real;
var
    s,s1 : Real;
begin
   if not key then // �������� ������� �� ����������
   begin
      my_count:=0.0;
      Laplas.MainMemo.Lines.Add('Error! expression can not be calculated... ');
      Application.MessageBox(PChar('expression can not be calculated...'),'Error!',MB_OK);
   end
   else
   begin
      if r^.myoperator[1]='o' then my_count:=r^.value // ���������
      else
       case r^.myoperator[1] of
         '+' :  begin
                    if (r^.left=nil) then
                   begin
                      // ������� ���� (���� ����������� ��� �� ����� �����������.
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
                      // ������� �����
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
                   if abs(s)<1e-10 then // ������������ ����
                   begin
                       my_count:=0.0;
                       key:=false;
                   end
                   else my_count:=my_count(r^.left,key)/s;
                end;
         '^' : begin
                  s:=my_count(r^.left,key);
                  s1:=my_count(r^.right,key);
                  if abs(s)<1e-10 then // ������������ ����
                  begin
                     if abs(s1)<1e-10 then my_count:=1.0
                     else  my_count:=0.0;
                  end
                  else
                  begin
                      my_count:=exp(s1*ln(abs(s)));
                  end;
               end;
         // ��������� �������� �������� �������������� �������
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
                     // sh - ��������������� �����
                     s:=my_count(r^.left,key);
                     my_count:=0.5*(exp(s)-exp(-s));
                  end;
                  if  r^.myoperator[2]='q' then
                  begin
                     if (length(r^.myoperator[2])=4) then
                     begin
                        // sqrt(abs(s)) ���������� ������
                        my_count:=sqrt(abs(my_count(r^.left,key)));
                     end
                     else
                     begin
                        // sqr - ���������� � �������
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
                     // ch ��������������� �������
                     s:=my_count(r^.left,key);
                     my_count:=0.5*(exp(s)+exp(-s));
                  end;
               end;
         'e' : begin
                   // exp   ����������
                   my_count:=exp(my_count(r^.left,key));
               end;
         'l' : begin
                  // ln(abs(s))  ����������� ��������
                  s:=my_count(r^.left,key);
                  my_count:=ln(abs(s));
               end;
         else // ������������� ��������
          begin
             my_count:=0.0;
             key:=false;
          end;
       end; // case

   end;
end; // my_count

// ������������ ����������� ������ ���������� �������� �������
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

// ����������� ���������������� ������ ���������� ���������� � ������������ �����.
function TFormVariables.my_real_convert(s : String; var bOk : Boolean) : Real;
var
   s1 : String;  // ������� ������
   Root : pTop; // ������ ������ ���������
   key : Boolean; // ������� ������������� ��������� � �������� �����
   r1 : Real; // ������������ ��������
   ileft, iright, i : Integer;

begin
    // ������� ����� ����������� � ���������� ������.
    ileft:=0;
    iright:=0;
    for i := 1 to length(s) do
    begin
       if (s[i]='(') then
       begin
           ileft:=ileft+1;
       end;
       if (s[i]=')') then
       begin
          iright:=iright+1;
       end;
    end;

    if (ileft=iright) then
    begin
       // ����� ����������� ������ ��������� ����� ����������� ������.


       // ����������� - ������ ���� ���������� �� ����������:
       my_substitutional_of_value(s ,s1);
       // ���� ���������� ��������, �� �� ���������� �������� � ������ ���������������� ����������.
       // ���������� ������.
       new(Root); // ��������� ����������� ������
       Constr_Tree(Root, s1); // ���������� ������ ���������
       key:=true;
       r1:=my_count(Root, key);
       my_delete(Root); // ������������ ����������� ������
       bOk:=key; // ������������ ������ ������� ���������� ��������.
       if not key then
       begin
          r1:=0.0;
          Laplas.MainMemo.Lines.Add('Error! expression :');
          Laplas.MainMemo.Lines.Add(s);
          Laplas.MainMemo.Lines.Add('can not be calculated... ');
          Application.MessageBox(PChar('expression '+s+' can not be calculated...'),'Error!',MB_OK);
       end;
       my_real_convert:=r1;
    end
     else
    begin
       // ���������� ����������� ������ �� ����� ����� ����������� ������.
       bOk:=false;
       my_real_convert:=-10000.0; // ���� ��������� ��������.
    end;
end;

// ���������� �������� ����� ��������� �������� ����������
// ������� ������ �������������.
procedure TFormVariables.my_update_size();
var
    i, i_4 : Integer;
    // ��������������� ���������� ��� ��������� �������������� ��������
    bOk : Boolean;
    r1, r2, r3, r4, r5, r6, rpow : Real;
    r7, r8, r9, r10, r11, r12 : Real;

begin
   // ������������� :
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

   // ����� :
   for i:=0 to (Laplas.lb-1) do
   begin
      // ���� �� ���� ������ � ��������� �������.
      bOk:=true; // ������� ������������ �����

      // ������������� �����������.

      if bOk then r1:=FormVariables.my_real_convert(Laplas.body[i].semissW,bOk);  // �������� �������
      if bOk then r2:=FormVariables.my_real_convert(Laplas.body[i].semissE,bOk);  // �������� �������������
      if bOk then r3:=FormVariables.my_real_convert(Laplas.body[i].semissS,bOk);  // � ������ �����������
      if bOk then r4:=FormVariables.my_real_convert(Laplas.body[i].semissN,bOk);  // �������� ����������.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.body[i].semissB,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.body[i].semissT,bOk);

      if (bOk)  then
      begin
        // ���������� �������� ������������� �����������.
       Laplas.body[i].emissW:=r1;
       Laplas.body[i].emissE:=r2;
       Laplas.body[i].emissS:=r3;
       Laplas.body[i].emissN:=r4;
       Laplas.body[i].emissB:=r5;
       Laplas.body[i].emissT:=r6;
      end;


      // ���������� �����

      if bOk then r1:=FormVariables.my_real_convert(Laplas.body[i].sxS,bOk);  // �������� �������
      if bOk then r2:=FormVariables.my_real_convert(Laplas.body[i].syS,bOk);  // �������� �������������
      if bOk then r3:=FormVariables.my_real_convert(Laplas.body[i].sxE,bOk);  // � ������ �����������
      if bOk then r4:=FormVariables.my_real_convert(Laplas.body[i].syE,bOk);  // �������� ����������.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.body[i].szS,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.body[i].szE,bOk);

        if bOk then r7:=FormVariables.my_real_convert(Laplas.body[i].sxC,bOk);  // �������� �������
      if bOk then r8:=FormVariables.my_real_convert(Laplas.body[i].syC,bOk);  // �������� �������������
      if bOk then r9:=FormVariables.my_real_convert(Laplas.body[i].szC,bOk);  // � ������ �����������
      if bOk then r10:=FormVariables.my_real_convert(Laplas.body[i].sHcyl,bOk);  // �������� ����������.
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
         Laplas.body[i].xS:=r1;  // �������� �������
         Laplas.body[i].yS:=r2;  // �������� �������������
         Laplas.body[i].xE:=r3;  // � ������ �����������
         Laplas.body[i].yE:=r4;  // �������� ����������.
         Laplas.body[i].zS:=r5;
         Laplas.body[i].zE:=r6;

         Laplas.body[i].xC:=r7;  // �������� �������
         Laplas.body[i].yC:=r8;  // �������� �������������
         Laplas.body[i].zC:=r9;  // � ������ �����������
         Laplas.body[i].Hcyl:=r10;  // �������� ����������.
         Laplas.body[i].R_out_cyl:=r11;
         Laplas.body[i].R_in_cyl:=r12;
         for i_4 := 0 to Laplas.body[i].n_power-1 do
         begin
             Laplas.body[i].arr_power[i_4]:=FormVariables.my_real_convert(Laplas.body[i].arr_s_power[i_4],bOk); // �������� ��������������
         end;
      end;


      Laplas.body[i].bactivity:=Laplas.activity_body(Laplas.body[i]);

   end;

   // ��������� �����
   for i:=0 to (Laplas.ls-1) do
   begin
      // ���� �� ���� ���������� ����� � ��������� �������.
      bOk:=true; // ������� ������������ �����

      // ���������� �����

      if bOk then r1:=FormVariables.my_real_convert(Laplas.source[i].sxS,bOk);  // �������� �������
      if bOk then r2:=FormVariables.my_real_convert(Laplas.source[i].syS,bOk);  // �������� �������������
      if bOk then r3:=FormVariables.my_real_convert(Laplas.source[i].sxE,bOk);  // � ������ �����������
      if bOk then r4:=FormVariables.my_real_convert(Laplas.source[i].syE,bOk);  // �������� ����������.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.source[i].szS,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.source[i].szE,bOk);
      if bOk then rpow:=FormVariables.my_real_convert(Laplas.source[i].spower,bOk);

      if (bOk) then
      begin
         Laplas.source[i].xS:=r1;  // �������� �������
         Laplas.source[i].yS:=r2;  // �������� �������������
         Laplas.source[i].xE:=r3;  // � ������ �����������
         Laplas.source[i].yE:=r4;  // �������� ����������.
         Laplas.source[i].zS:=r5;
         Laplas.source[i].zE:=r6;
         Laplas.source[i].Power:=rpow; // �������� ��������������

         Laplas.source[i].bactivity:=Laplas.activity_plane(Laplas.source[i]);
      end;
   end;

   // ������ ������
   for i:=0 to (Laplas.lw-1) do
   begin
      // ���� �� ���� ���������� ����� � ��������� �������.
      bOk:=true; // ������� ������������ �����

      // ���������� �����

      if bOk then r1:=FormVariables.my_real_convert(Laplas.wall[i].sxS,bOk);  // �������� �������
      if bOk then r2:=FormVariables.my_real_convert(Laplas.wall[i].syS,bOk);  // �������� �������������
      if bOk then r3:=FormVariables.my_real_convert(Laplas.wall[i].sxE,bOk);  // � ������ �����������
      if bOk then r4:=FormVariables.my_real_convert(Laplas.wall[i].syE,bOk);  // �������� ����������.
      if bOk then r5:=FormVariables.my_real_convert(Laplas.wall[i].szS,bOk);
      if bOk then r6:=FormVariables.my_real_convert(Laplas.wall[i].szE,bOk);

      if (bOk) then
      begin
         Laplas.wall[i].xS:=r1;  // �������� �������
         Laplas.wall[i].yS:=r2;  // �������� �������������
         Laplas.wall[i].xE:=r3;  // � ������ �����������
         Laplas.wall[i].yE:=r4;  // �������� ����������.
         Laplas.wall[i].zS:=r5;
         Laplas.wall[i].zE:=r6;
      end;

       bOk:=true; // ������� ������������ �����
       if bOk then r1:=FormVariables.my_real_convert(Laplas.wall[i].semissivity,bOk);  // �������� �������
       if bOk then r2:=FormVariables.my_real_convert(Laplas.wall[i].sheat_transfer_coefficient,bOk);

      if (bOk) then
      begin
         // ������������� �����������.
         Laplas.wall[i].emissivity:=r1;
         // ����������� �����������.
         Laplas.wall[i].heat_transfer_coefficient:=r2;

         Laplas.wall[i].bactivity:=Laplas.activity_plane(Laplas.wall[i]);
      end;

   end;

end;

end.
