unit UnitOrientation_pcb_plate_in_space;
// pcb ����� ������ ����� � ��������� XOY �� ������ ����� ���� ����������
// ��� � ������� � ������� ������ OZ ��� � � �������� OY.
// � ������ ������� OY ������������� �������������� ����� XOY->XOZ.
// 29.11.2022

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormOrientation_pcb_plate = class(TForm)
    RadioGroupOrientation_pcb: TRadioGroup;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormOrientation_pcb_plate: TFormOrientation_pcb_plate;

implementation

{$R *.dfm}

end.
