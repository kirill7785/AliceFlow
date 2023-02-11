unit UnitFormUnion;
// ���������� ��������

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormUnion = class(TForm)
    PanelMain: TPanel;
    ButtonApply: TButton;
    CheckBoxVisible: TCheckBox;
    Labelname: TLabel;
    Editname: TEdit;
    CBMeshassemblesSeparat: TCheckBox;
    GBmeshasembles: TGroupBox;
    GBMesh: TGroupBox;
    Linx: TLabel;
    CBinx: TComboBox;
    Liny: TLabel;
    Linz: TLabel;
    CBiny: TComboBox;
    CBinz: TComboBox;
    GBadditsize: TGroupBox;
    LXmin: TLabel;
    LYmin: TLabel;
    LZmin: TLabel;
    Exmin: TEdit;
    Eymin: TEdit;
    Ezmin: TEdit;
    LXmax: TLabel;
    LYmax: TLabel;
    LZmax: TLabel;
    EXmax: TEdit;
    EYmax: TEdit;
    EZmax: TEdit;
    Labelmassa: TLabel;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    ComboBoxUnionMeshGenerator: TComboBox;
    procedure ButtonApplyClick(Sender: TObject);
    procedure CheckBoxVisibleClick(Sender: TObject);
    procedure CBMeshassemblesSeparatClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormUnion: TFormUnion;

implementation
 uses
     VisualUnit, UnitVariables;

{$R *.dfm}


// ��������� �������� � �����������
procedure TFormUnion.ButtonApplyClick(Sender: TObject);
var
    bvisibleunite  : Boolean;
    i  : Integer;
    bOk : Boolean;
    rbuf : Real;
begin
   bvisibleunite:=CheckBoxVisible.Checked;
   CheckBoxVisibleClick(Sender);

   Laplas.myassembles[Laplas.itek-1].name:=Editname.Text; // ����� ��� �����������.
   // ������������� ����� ������� ����� �������� ����������� ���.
   Laplas.correctobjname('u',Laplas.myassembles[Laplas.itek-1].name,Laplas.itek-1);
   Editname.Text:=Laplas.myassembles[Laplas.itek-1].name;
   FormUnion.Caption:='Assemblies [ '+Trim(Laplas.myassembles[Laplas.itek-1].name)+' ]';

   // ��������� ������ �������� � ����������� �������������� ��� ����� �� ������ ��������� �� �����.

   Laplas.myassembles[Laplas.itek-1].bmesh_assembles_separately:=CBMeshassemblesSeparat.Checked;
   if (Laplas.myassembles[Laplas.itek-1].bmesh_assembles_separately) then
   begin
      bOk:=true;
      if bOk then rbuf:=FormVariables.my_real_convert(Trim(Exmin.Text),bOk);
      if bOk then rbuf:=FormVariables.my_real_convert(Trim(Eymin.Text),bOk);
      if bOk then rbuf:=FormVariables.my_real_convert(Trim(Ezmin.Text),bOk);
      if bOk then rbuf:=FormVariables.my_real_convert(Trim(EXmax.Text),bOk);
      if bOk then rbuf:=FormVariables.my_real_convert(Trim(EYmax.Text),bOk);
      if bOk then rbuf:=FormVariables.my_real_convert(Trim(EZmax.Text),bOk);
      if (bOk) then
      begin
         // ����������������� ������� �����.
         Laplas.myassembles[Laplas.itek-1].sxmin:=Trim(Exmin.Text);
         Laplas.myassembles[Laplas.itek-1].symin:=Trim(Eymin.Text);
         Laplas.myassembles[Laplas.itek-1].szmin:=Trim(Ezmin.Text);
         Laplas.myassembles[Laplas.itek-1].sxmax:=Trim(EXmax.Text);
         Laplas.myassembles[Laplas.itek-1].symax:=Trim(EYmax.Text);
         Laplas.myassembles[Laplas.itek-1].szmax:=Trim(EZmax.Text);
         // �������� ������� �����.
         Laplas.myassembles[Laplas.itek-1].xmin:=FormVariables.my_real_convert(Trim(Exmin.Text),bOk);
         Laplas.myassembles[Laplas.itek-1].ymin:=FormVariables.my_real_convert(Trim(Eymin.Text),bOk);
         Laplas.myassembles[Laplas.itek-1].zmin:=FormVariables.my_real_convert(Trim(Ezmin.Text),bOk);
         Laplas.myassembles[Laplas.itek-1].xmax:=FormVariables.my_real_convert(Trim(EXmax.Text),bOk);
         Laplas.myassembles[Laplas.itek-1].ymax:=FormVariables.my_real_convert(Trim(EYmax.Text),bOk);
         Laplas.myassembles[Laplas.itek-1].zmax:=FormVariables.my_real_convert(Trim(EZmax.Text),bOk);
         //  �������� ���������.
         Laplas.myassembles[Laplas.itek-1].inxloc:=CBinx.ItemIndex+4;
         Laplas.myassembles[Laplas.itek-1].inyloc:=CBiny.ItemIndex+4;
         Laplas.myassembles[Laplas.itek-1].inzloc:=CBinz.ItemIndex+4;
         Laplas.myassembles[Laplas.itek-1].itypeMesh:=ComboBoxUnionMeshGenerator.ItemIndex;
      end
       else
      begin
         ShowMessage('Input incorrect. Please repeat once more...');
      end;
   end;
end;

procedure TFormUnion.CheckBoxVisibleClick(Sender: TObject);
var
    bvisibleunite  : Boolean;
    i, j, ident_now  : Integer;
    stack : array of Integer;
    top : Integer;

procedure my_visible_check(identifire : Integer; bvisibleunite_now : Boolean);
var
  i_1 : Integer;
begin
   // �����
   for i_1:=1 to (Laplas.lb-1) do
   begin
      if (Laplas.body[i_1].iunion=identifire) then Laplas.body[i_1].bvisible:=bvisibleunite_now;
   end;
   // ���������
   for i_1:=0 to (Laplas.ls-1) do
   begin
      if (Laplas.source[i_1].iunion=identifire) then Laplas.source[i_1].bvisible:=bvisibleunite_now;
   end;
   // ������
   for i_1:=0 to (Laplas.lw-1) do
   begin
      if (Laplas.wall[i_1].iunion=identifire) then Laplas.wall[i_1].bvisible:=bvisibleunite_now;
   end;
end;


begin
   // ��������� ��� ������ ������ ���������.

   top:=0;
   SetLength(stack,Laplas.lu);
   stack[top]:= Laplas.myassembles[Laplas.itek-1].identifire;

   bvisibleunite:=CheckBoxVisible.Checked;
   Laplas.myassembles[Laplas.itek-1].bVisible:=bvisibleunite;
   my_visible_check(Laplas.myassembles[Laplas.itek-1].identifire, bvisibleunite);

   while (top>-1) do
   begin

      ident_now:=stack[top];
      dec(top);

      for j := 0 to Laplas.lu-1 do
      begin
         if (Laplas.myassembles[j].iunionparent+1=ident_now) then
         begin
            inc(top);
            stack[top]:=Laplas.myassembles[j].identifire;

            Laplas.myassembles[j].bVisible:=CheckBoxVisible.Checked;

            my_visible_check(Laplas.myassembles[j].identifire, bvisibleunite);
         end;
      end;



   end;

   Laplas.ReadyPaint;
end;

// ��������� ��� ���������� ������ ����������������� �����.
procedure TFormUnion.CBMeshassemblesSeparatClick(Sender: TObject);
begin
   Laplas.myassembles[Laplas.itek-1].bmesh_assembles_separately:=CBMeshassemblesSeparat.Checked;
   GBmeshasembles.Visible:=CBMeshassemblesSeparat.Checked;
   Exmin.Text:=Laplas.myassembles[Laplas.itek-1].sxmin;
   Eymin.Text:=Laplas.myassembles[Laplas.itek-1].symin;
   Ezmin.Text:=Laplas.myassembles[Laplas.itek-1].szmin;
   EXmax.Text:=Laplas.myassembles[Laplas.itek-1].sxmax;
   EYmax.Text:=Laplas.myassembles[Laplas.itek-1].symax;
   EZmax.Text:=Laplas.myassembles[Laplas.itek-1].szmax;
   // ����� :
   CBinx.ItemIndex:=Laplas.myassembles[Laplas.itek-1].inxloc-4;
   CBiny.ItemIndex:=Laplas.myassembles[Laplas.itek-1].inyloc-4;
   CBinz.ItemIndex:=Laplas.myassembles[Laplas.itek-1].inzloc-4;
end;

end.
