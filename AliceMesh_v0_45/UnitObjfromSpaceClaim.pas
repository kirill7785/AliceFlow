unit UnitObjfromSpaceClaim;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormObjfromSpaceClaim = class(TForm)
    ButtonApply: TButton;
    GroupBoxextrudesource2D23D: TGroupBox;
    LabelX: TLabel;
    Label1: TLabel;
    Label2: TLabel;
    ComboBoxnormalX: TComboBox;
    ComboBoxnormalY: TComboBox;
    ComboBoxnormalZ: TComboBox;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormObjfromSpaceClaim: TFormObjfromSpaceClaim;

implementation

{$R *.dfm}

procedure TFormObjfromSpaceClaim.ButtonApplyClick(Sender: TObject);
begin
   Close;
end;

end.
