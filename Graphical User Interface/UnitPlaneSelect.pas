unit UnitPlaneSelect;
// Выбор плоскости при повороте.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormSelectPlaneRotation = class(TForm)
    Label1: TLabel;
    ComboBoxPlane: TComboBox;
    ButtonApply: TButton;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormSelectPlaneRotation: TFormSelectPlaneRotation;

implementation

{$R *.dfm}

procedure TFormSelectPlaneRotation.ButtonApplyClick(Sender: TObject);
begin
    Close;
end;

end.
