unit UnitTextNameSourcePattern;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormTextNameSourcePattern = class(TForm)
    Label1: TLabel;
    EditSourcepatternname: TEdit;
    ButtonApply: TButton;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormTextNameSourcePattern: TFormTextNameSourcePattern;

implementation

{$R *.dfm}

procedure TFormTextNameSourcePattern.ButtonApplyClick(Sender: TObject);
begin
   EditSourcepatternname.Text:=Trim(EditSourcepatternname.Text);
   if (length(EditSourcepatternname.Text)>0) then
   begin
      // Переходим к заданию тепловой мощности в образцах.
      EditSourcepatternname.Color:=clWhite;
      Close;
   end
   else
   begin
      EditSourcepatternname.Color:=clRed;
   end;
end;

end.
