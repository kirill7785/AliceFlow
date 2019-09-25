unit UnitAMGCLManager;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormAMGCLParameters = class(TForm)
    RadioGroupAMGCLsmoother1: TRadioGroup;
    Button1: TButton;
    RadioGroupAMGCLCoarseningType: TRadioGroup;
    Label1: TLabel;
    ComboBoxIterator: TComboBox;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormAMGCLParameters: TFormAMGCLParameters;

implementation

{$R *.dfm}

procedure TFormAMGCLParameters.Button1Click(Sender: TObject);
begin
   Close;
end;

end.
