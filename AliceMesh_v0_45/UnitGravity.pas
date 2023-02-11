unit UnitGravity;
// 6 августа 2016 данный модуль устарел и его функции
// теперь исполняет UnitEQGD Form.
// Т.к. UnitEQGD Form теперь имеет вид как форма Basic Parameters
// в ANSYS Icepak 12.0, что позволило упростить логику работы пользователя
// с программой AliceMesh_v0_38 && AliceFlow_v0_24.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls;

type
  TFormGravity = class(TForm)
    GBGravity: TGroupBox;
    Lgx: TLabel;
    Lgy: TLabel;
    Lgz: TLabel;
    Egx: TEdit;
    Egy: TEdit;
    Egz: TEdit;
    Bapply: TButton;
    lblgx: TLabel;
    lblgy: TLabel;
    lbldz: TLabel;
    procedure BapplyClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
    procedure patchstring(var s : String);
  public
    { Public declarations }
  end;

var
  FormGravity: TFormGravity;

implementation
 uses
     VisualUnit;
{$R *.dfm}


procedure TFormGravity.patchstring(var s : String);
var
   i : Integer;
begin
if (FormatSettings.DecimalSeparator='.') then
begin
   for i:=1 to length(s) do
   begin
      if (s[i]=',') then s[i]:='.';
   end;
end;

if (FormatSettings.DecimalSeparator=',') then
begin
   for i:=1 to length(s) do
   begin
      if (s[i]='.') then s[i]:=',';
   end;
end;
end;

// Ввод значения силы тяжести по нажатию на кнопку Apply-применить.
procedure TFormGravity.BapplyClick(Sender: TObject);
var
  s : String;
  bOk : Boolean;
  c : Real;
  code : Integer;
  sforval : String;
begin

   // исправление десятичного разделителя.
   s:=Trim(Egx.Text);
   patchstring(s);
   Egx.Text:=s;
   s:=Trim(Egy.Text);
   patchstring(s);
   Egy.Text:=s;
   s:=Trim(Egz.Text);
   patchstring(s);
   Egz.Text:=s;


   bOk:=False;
   //val(Trim(Egx.Text),c,code);
   sforval:='';
  sforval:=StringReplace(Egx.Text,',','.',[rfReplaceAll]);
   val(sforval,c,code);
   if (code=0) then
   begin
       //val(Trim(Egy.Text),c,code);
       sforval:='';
       sforval:=StringReplace(Egy.Text,',','.',[rfReplaceAll]);
       val(sforval,c,code);
      if (code=0) then
      begin
         //val(Trim(Egz.Text),c,code);
         sforval:='';
         sforval:=StringReplace(Egz.Text,',','.',[rfReplaceAll]);
         val(sforval,c,code);
         if (code=0) then
         begin
            bOk:=True;
         end;
      end;
   end;
   if (bOk) then
   begin
      Laplas.gx:=StrToFloat(Trim(Egx.Text));
      Laplas.gy:=StrToFloat(Trim(Egy.Text));
      Laplas.gz:=StrToFloat(Trim(Egz.Text));
   end
    else
   begin
      if (FormatSettings.DecimalSeparator='.') then
      begin
         Egx.Text:='0.0';
         Egy.Text:='0.0';
         Egz.Text:='0.0';
      end;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         Egx.Text:='0,0';
         Egy.Text:='0,0';
         Egz.Text:='0,0';
      end;
      ShowMessage('Input Error, please repeat...');
   end;
end;

// Инициализация формы GravityForm при создании формы.
procedure TFormGravity.FormCreate(Sender: TObject);
begin
    Egx.Text:=FloatToStr(Laplas.gx);
    Egy.Text:=FloatToStr(Laplas.gy);
    Egz.Text:=FloatToStr(Laplas.gz);
end;

end.
