unit UnitExportSYMMIC;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls;

type
  TFormexport2SYMMIC = class(TForm)
    lbl1: TLabel;
    edtname: TEdit;
    lbl2: TLabel;
    edttitle: TEdit;
    lblnx: TLabel;
    lblny: TLabel;
    lblnz: TLabel;
    cbbnx: TComboBox;
    cbbny: TComboBox;
    cbbnz: TComboBox;
    lbltol: TLabel;
    cbbtol: TComboBox;
    btnApply: TButton;
    procedure btnApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Formexport2SYMMIC: TFormexport2SYMMIC;

implementation

{$R *.dfm}

procedure TFormexport2SYMMIC.btnApplyClick(Sender: TObject);
var
  s : String;
  i : Integer;
  flag : Boolean;
begin
   // убираем лишние пробелы.
   edtname.Text:=Trim(edtname.Text);
   edttitle.Text:=Trim(edttitle.Text);
   if ((Length(edtname.Text)>0) and (Length(edttitle.Text)>0)) then
   begin
      flag:=true;
      s:='';
      for i:=1 to length(edtname.Text) do
      begin
          if (edtname.Text[i]='.') then
          begin
             flag:=false;
          end;
          if (flag) then
          begin
             s[i]:=edtname.Text[i];
          end;
      end;
      edtname.Text:=s; // обрезали всё что после точки.
      if ((Length(edtname.Text)>0) and (Length(edttitle.Text)>0)) then
      begin
         Close;
      end;
   end
   else
   begin
      // не введены строки имени файла и заголовка.
      ShowMessage('Invalid string file name and title.');
   end;
end;



end.
