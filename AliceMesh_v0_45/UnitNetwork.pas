unit UnitNetwork;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormNetwork = class(TForm)
    GroupBox1: TGroupBox;
    CheckBoxNetworkPartition: TCheckBox;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormNetwork: TFormNetwork;

implementation

{$R *.dfm}

end.
