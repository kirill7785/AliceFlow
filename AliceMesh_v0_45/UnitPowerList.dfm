object FormPowerList: TFormPowerList
  Left = 322
  Top = 114
  AutoSize = True
  Caption = 'specify list table'
  ClientHeight = 321
  ClientWidth = 209
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Panelmain: TPanel
    Left = 0
    Top = 0
    Width = 209
    Height = 321
    Color = clMoneyGreen
    TabOrder = 0
    object Label1: TLabel
      Left = 16
      Top = 16
      Width = 142
      Height = 13
      Caption = 'Please, specify a list tabulated'
    end
    object Label2: TLabel
      Left = 16
      Top = 40
      Width = 132
      Height = 13
      Caption = 'by Pdiss as a function Tmax'
    end
    object Label3: TLabel
      Left = 16
      Top = 64
      Width = 76
      Height = 13
      Caption = 'and offset drain.'
    end
    object Label4: TLabel
      Left = 16
      Top = 88
      Width = 78
      Height = 13
      Caption = 'number of tables'
    end
    object Eltdp: TEdit
      Left = 104
      Top = 80
      Width = 65
      Height = 21
      TabOrder = 0
    end
    object BApplyglobal: TButton
      Left = 96
      Top = 112
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = BApplyglobalClick
    end
    object GBtableparam: TGroupBox
      Left = 16
      Top = 144
      Width = 177
      Height = 161
      Caption = 'Table parameters'
      TabOrder = 2
      object Label5: TLabel
        Left = 16
        Top = 32
        Width = 34
        Height = 13
        Caption = 'id table'
      end
      object Label6: TLabel
        Left = 16
        Top = 64
        Width = 68
        Height = 13
        Caption = 'table file name'
      end
      object CBidtable: TComboBox
        Left = 64
        Top = 24
        Width = 57
        Height = 21
        TabOrder = 0
        OnChange = CBidtableChange
      end
      object Efilename: TEdit
        Left = 16
        Top = 88
        Width = 145
        Height = 21
        TabOrder = 1
      end
      object Bonetable: TButton
        Left = 80
        Top = 120
        Width = 75
        Height = 25
        Caption = 'Apply'
        TabOrder = 2
        OnClick = BonetableClick
      end
    end
  end
end
