unit UnitCopy;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormCopyObject = class(TForm)
    PanelMain: TPanel;
    Labelnc: TLabel;
    Editnc: TEdit;
    GroupBoxTranslate: TGroupBox;
    Labelx: TLabel;
    Labely: TLabel;
    Labelz: TLabel;
    EditX: TEdit;
    EditY: TEdit;
    EditZ: TEdit;
    ButtonApply: TButton;
    CheckBoxRotate: TCheckBox;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    procedure ButtonApplyClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormClose(Sender: TObject; var Action: TCloseAction);
  private
    { Private declarations }
    bvisit : Boolean;
    procedure patchstring(var s : String);
  public
    { Public declarations }
  end;

var
  FormCopyObject: TFormCopyObject;

implementation
uses
     VisualUnit, UnitVariables;
{$R *.dfm}

procedure TFormCopyObject.patchstring(var s : String);
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

procedure TFormCopyObject.ButtonApplyClick(Sender: TObject);
var
    bOk : Boolean;
    s : String;
begin
   // ввод значений в основную форму
   Laplas.inum:=StrToInt(Editnc.Text); // число копий

   s:=Trim(EditX.Text);
   patchstring(s);
   EditX.Text:=s;
   s:=Trim(EditY.Text);
   patchstring(s);
   EditY.Text:=s;
   s:=Trim(EditZ.Text);
   patchstring(s);
   EditZ.Text:=s;

   // проверка корректности ввода.
   bOk:=true;
   FormVariables.my_real_convert(EditX.Text,bOk);
   if (bOk) then FormVariables.my_real_convert(EditY.Text,bOk);
   if (bOk) then FormVariables.my_real_convert(EditZ.Text,bOk);



   Laplas.bcontinuecopy:=bOk;

   if (bOk) then
   begin
      Laplas.spdx:=EditX.Text; // параметризованное смещение по оси X
      Laplas.spdy:=EditY.Text; // параметризованное смещение по оси Y
      Laplas.spdz:=EditZ.Text; // параметризованное смещение по оси Z
      
      bvisit:=true;
      Close;
   end;
end;



procedure TFormCopyObject.FormCreate(Sender: TObject);
begin
   bvisit:=false;
end;

procedure TFormCopyObject.FormClose(Sender: TObject;
  var Action: TCloseAction);
begin
   if (bvisit) then
   begin
      bvisit:=false;
      Close;
   end
    else
   begin
       Laplas.inum:=0;
       Laplas.bcontinuecopy:=false;
       Close;
   end;
end;

end.
