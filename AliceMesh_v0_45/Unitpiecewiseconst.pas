unit Unitpiecewiseconst;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormpiecewiseconstant = class(TForm)
    Memopiecewiseconst: TMemo;
    Label2: TLabel;
    Label3: TLabel;
    ComboBoxpiecewiseconst: TComboBox;
    ButtonClear: TButton;
    procedure ButtonClearClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Formpiecewiseconstant: TFormpiecewiseconstant;

implementation

{$R *.dfm}

// ќчистка меню редактировани€ piecewise constant шагов по времени.
procedure TFormpiecewiseconstant.ButtonClearClick(Sender: TObject);
begin
   Memopiecewiseconst.Clear;
end;

end.
