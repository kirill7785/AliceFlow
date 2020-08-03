unit UnitScale;
// 04.09.2019 Модуль в котором можно проивести
// масштабирование всей пользовательской геометрии.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormScale = class(TForm)
    GroupBox1: TGroupBox;
    Label1: TLabel;
    Edit1: TEdit;
    ButtonScale: TButton;
    procedure ButtonScaleClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    scale674 : Real;
  end;

var
  FormScale: TFormScale;

implementation

{$R *.dfm}

procedure TFormScale.ButtonScaleClick(Sender: TObject);
var
  s : String;
  code : Integer;
  dmbuf : Real;

begin
   s:=Edit1.Text;
   if (FormatSettings.DecimalSeparator=',') then
   begin
      s:=StringReplace(s,'.',',',[rfReplaceAll]);
   end;
   if (FormatSettings.DecimalSeparator='.') then
   begin
      s:=StringReplace(s,',','.',[rfReplaceAll]);
   end;
   val(s,dmbuf,code);
   if (code=0) then
   begin
      if (abs(dmbuf)<1.0e-30) then
      begin
          // Защита от нуля.
          dmbuf:=1.0;
      end;
      if (dmbuf>1.0e-30) then
      begin
         scale674:=abs(dmbuf); // защита от отрицательного значения.
         Edit1.Color:=clWhite;
         // Автоматическое закрытие формы.
         Close;
      end
      else
      begin
         Edit1.Color:=clRed;
      end;
   end
    else
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         Edit1.Text:='1,0';
      end
      else
      begin
         Edit1.Text:='1.0';
      end;
   end;
end;

procedure TFormScale.FormCreate(Sender: TObject);
begin
   scale674:=1.0; // Инициализация при создании формы 19.05.2020.
end;

end.
