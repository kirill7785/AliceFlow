unit Unitamg1r5Parameters;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormamg1r5Parameters = class(TForm)
    GroupBoxStabilisation: TGroupBox;
    ComboBoxStabilization: TComboBox;
    ButtonApply: TButton;
    CheckBox_amg1r6: TCheckBox;
    GroupBox1: TGroupBox;
    ComboBoxNumber_of_smootherssteps: TComboBox;
    ComboBoxTypeSmoother: TComboBox;
    Labelpre: TLabel;
    LabelpostSmooth: TLabel;
    ComboBoxNumber_of_post_smooth: TComboBox;
    ComboBoxTypePostSmoother: TComboBox;
    RGthresholds: TRadioGroup;
    Labelstrongthreshold: TLabel;
    Editstrongthreshold: TEdit;
    LabelF2F: TLabel;
    EditF2F: TEdit;
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
