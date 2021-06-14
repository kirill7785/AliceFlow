unit UnitOptimetric;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormOptimetric = class(TForm)
    EditListVariable: TEdit;
    Label1: TLabel;
    ComboBoxvar_id0: TComboBox;
    ButtonRun: TButton;
    ComboBoxvar_id1: TComboBox;
    EditListVariable1: TEdit;
    procedure ButtonRunClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }

  end;

var
  FormOptimetric: TFormOptimetric;

implementation

{$R *.dfm}

uses VisualUnit;

procedure TFormOptimetric.ButtonRunClick(Sender: TObject);
begin
   Laplas.bOkTrials:=true;
   Close;
end;

end.
