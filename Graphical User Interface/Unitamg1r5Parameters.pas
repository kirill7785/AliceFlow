unit Unitamg1r5Parameters;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormamg1r5Parameters = class(TForm)
    GroupBoxStabilisation: TGroupBox;
    ComboBoxStabilization: TComboBox;
    ButtonApply: TButton;
    CheckBox_amg1r6: TCheckBox;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Formamg1r5Parameters: TFormamg1r5Parameters;

implementation

{$R *.dfm}

procedure TFormamg1r5Parameters.ButtonApplyClick(Sender: TObject);
begin
   Close;
end;

end.
