unit UnitPatternDelete;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormPattern = class(TForm)
    RadioGroupPattern: TRadioGroup;
    Label1: TLabel;
    EdittextnamefragmentPattern: TEdit;
    ButtonApply: TButton;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormPattern: TFormPattern;

implementation

{$R *.dfm}

procedure TFormPattern.ButtonApplyClick(Sender: TObject);
begin
   EdittextnamefragmentPattern.Text:=Trim(EdittextnamefragmentPattern.Text);
   Close;
end;

end.
