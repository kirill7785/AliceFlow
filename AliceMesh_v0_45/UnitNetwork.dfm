object FormNetwork: TFormNetwork
  Left = 0
  Top = 0
  Caption = 'Network partition'
  ClientHeight = 242
  ClientWidth = 422
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 8
    Top = 8
    Width = 169
    Height = 57
    Caption = 'partition'
    TabOrder = 0
    object CheckBoxNetworkPartition: TCheckBox
      Left = 32
      Top = 24
      Width = 97
      Height = 17
      Caption = 'Network partition'
      Checked = True
      State = cbChecked
      TabOrder = 0
    end
  end
end
